struct bcol_model_t {
  fan_3d::model::fms_t fms;
  BCOL_t::ObjectID_t oid;

  bcol_model_t() = default;
  bcol_model_t(const char* name) : fms({ .path = name }) {

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
    g_bcol.ClearObject(oid);
    static f32_t totald = 0;
    
    totald += d * 150;
    fms.dt = totald;
    //fms.calculate_default_pose();
    //fms.fk_get_pose(totald);

    //auto fk_animation_transform = fms.calculate_transformations();
    //hardcode
    //static constexpr int mesh_id = 0;
    //for (uintptr_t i = 0; i < fms.parsed_model.model_data.mesh_data.size(); i++) {
      //fms.calculate_modified_vertices(i, fk_animation_transform);
    //}

    struct triangle_list_t {
      //uint32_t matid; todo
      std::vector<fan_3d::model::fms_t::one_triangle_t> triangle_vec;
    };

    std::vector<triangle_list_t> triangles;
    for (uintptr_t i = 0; i < fms.parsed_model.model_data.mesh_data.size(); i++) {
      triangle_list_t tl;
      //tl.matid = fms.get_material_id(i);
      fms.get_triangle_vec(i, &tl.triangle_vec);

      triangles.push_back(tl);
    }
    for (uintptr_t i = 0; i < triangles.size(); i++) {
      auto& tl = triangles[i];
      for (uintptr_t ti = 0; ti < tl.triangle_vec.size(); ti++) {
        auto& t = tl.triangle_vec[ti];
        BCOL_t::ShapeProperties_DPF_t sp;
        sp.u.bcol_model = this;
        auto found = fan_3d::model::cached_texture_data.find(fms.parsed_model.model_data.mesh_data[i].names[aiTextureType_DIFFUSE]);
        if (found == fan_3d::model::cached_texture_data.end()) {
          sp.u.MaterialIndex = 0;
        }
        else {
          sp.u.MaterialIndex = std::distance(
            fan_3d::model::cached_texture_data.begin(),
            found
          ); // TODO
        }

        sp.u.color = fan::vec4(1); // TODO
        for (uint8_t pi = 0; pi < 3; pi++) {
          auto v = fms.m_transform * fan::vec4(t.p[pi], 1);
          sp.p[pi] = *(fan::vec3*)&v;
          sp.u.uv[pi] = t.tc[pi];
        }
        g_bcol.NewShape_DPF(oid, &sp);
      }
    }
  }
};
