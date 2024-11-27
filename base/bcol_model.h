struct bcol_model_t : fan_3d::model::fms_t{
  BCOL_t::ObjectID_t oid;

  /*for gui*/
  fan::vec3 model_position = 0;
  fan::vec3 model_rotation = fan::vec3(0, 0, 0);
  fan::vec3 model_scale = 0.1;

  bcol_model_t() = default;
  bcol_model_t(const char* name) : fms_t({ .path = name, .use_cpu = 1 }) {
    for (const fan_3d::model::mesh_t& mesh : meshes) {
      model_material_t material;
      auto& ci = fan_3d::model::cached_texture_data;
      auto found = ci.find(mesh.texture_names[aiTextureType_DIFFUSE]);
      if (found == ci.end()) {
        material.diffuse_id = 0;
      }
      else {
        material.diffuse_id = std::distance(ci.begin(), found);
      }
      model_material.push_back(material);
    }

  }
  void open() {
    BCOL_t::ObjectProperties_t p;
    //p.Position = BCOL_t::_vf(8, -115, 8);
    p.Position = BCOL_t::_vf(0, 0, 0);
    oid = g_bcol.NewObject(&p, BCOL_t::ObjectFlag::Constant);
  }
  void close() {
    if (oid.iic()) {
      return;
    }
    g_bcol.UnlinkObject(oid);
    g_bcol.RecycleObject(oid);
    oid.sic();
  }

  void delta(f32_t d) {
    ImGui::Begin("window");
    ImGui::DragFloat3("translate", model_position.data(), 0.1);
    ImGui::DragFloat3("rotate", model_rotation.data(), 0.01);
    ImGui::DragFloat3("scale", model_scale.data(), 0.01);

    g_bcol.ClearObject(oid);
    static f32_t totald = 0;
    
    totald += d * 1000;
    dt = totald;

    fk_calculate_poses();
    std::vector<fan::mat4> fk_transformations = bone_transforms;
    fk_interpolate_animations(fk_transformations, *root_bone, m_transform);

    mouse_modify_joint();

    fan::mat4 model_transform{1};
    model_transform = model_transform.translate(model_position);
    model_transform = model_transform.rotate(model_rotation);
    model_transform = model_transform.scale(model_scale);

    for (uint32_t i = 0; i < meshes.size(); ++i) {
      BCOL_t::ShapeProperties_DPF_t sp;
      sp.u.MaterialIndex = i;
      sp.u.bcol_model = this;
      calculate_vertices(fk_transformations, i, fan::mat4(1));

      for (size_t j = 0; j < meshes[i].indices.size(); j += 3) {
        for (uint8_t pi = 0; pi < 3; pi++) {
          uint32_t vertex_index = meshes[i].indices[j + pi];
          fan::vec4 v(model_transform * fan::vec4(calculated_meshes[i].vertices[vertex_index].position, 1.0));

          sp.p[pi] = *(fan::vec3*)&v;
          sp.u.uv[pi] = calculated_meshes[i].vertices[vertex_index].uv;
        }
        g_bcol.NewShape_DPF(oid, &sp);
      }
    }
    ImGui::End();
  }
};
