#ifndef RenderSizeX
  #define RenderSizeX 512
#endif
#ifndef RenderSizeY
  #define RenderSizeY 512
#endif
#ifndef RenderSizeDivide
  #define RenderSizeDivide 16
#endif
#ifndef set_PrintStats
  #define set_PrintStats 0
#endif

#ifndef set_VisualDebug
  #define set_VisualDebug 1
#endif
#ifndef set_HaveGrid
  #define set_HaveGrid 0
#endif
#ifndef set_DisplayWindow
  #define set_DisplayWindow 1
#endif
#ifndef set_stdout_data
  #define set_stdout_data 0
#endif
#ifndef set_RenderFrameCount
  #define set_RenderFrameCount (100000)
#endif
#ifndef set_RenderTick
  #define set_RenderTick 0.04
#endif

#ifndef set_bcol_UseEmbree
  #define set_bcol_UseEmbree 1
#endif
#ifndef set_bcol_IterateRaysMulti
  #define set_bcol_IterateRaysMulti 0
#endif
#ifndef set_Multithread
  #define set_Multithread 1
#endif
#ifndef set_Multithread_UseCond
  #define set_Multithread_UseCond 1
#endif

#if set_bcol_UseEmbree
  #include <embree4/rtcore.h>
#endif

#include <WITCH/WITCH.h>
#include <WITCH/IO/IO.h>
#include <WITCH/IO/print.h>
#include <WITCH/A/A.h>
#include <WITCH/T/T.h>

#include <fan/pch.h>
#include <fan/graphics/opengl/3D/objects/model.h>

void print(const char* format, ...) {
  IO_fd_t fd_stdout;
  IO_fd_set(&fd_stdout, FD_OUT);
  va_list argv;
  va_start(argv, format);
  IO_vprint(&fd_stdout, format, argv);
  va_end(argv);
}

template <typename ...Args>
constexpr void fprint(const Args&... args) {
  fan::print(args...);
}

#if set_Multithread == 1
  #include <WITCH/TH/TH.h>
#endif

fan::vec2ui RenderSize = { RenderSizeX, RenderSizeY };

constexpr f32_t GridBlockSize = 32;

bool IsGridSolid(fan::vec3si gi) {
  if (gi[1] == -5 && gi[0] == 0 && gi[2] == 0) {
    return true;
  }
  for (uint32_t grid_d = 0; grid_d < 3; grid_d++) {
    if (abs(gi[grid_d]) > 5) {
      return true;
    }
  }
  return false;
}

#if set_VisualDebug == 1
  void DebugBall(fan::vec3 p);
#endif

struct bcol_model_t;

#define BCOL_set_VisualSolve_dmin 0.001
#define BCOL_set_prefix BCOL
#define BCOL_set_SupportGrid set_HaveGrid
#define BCOL_set_Dimension 3
#define BCOL_set_StoreExtraDataInsideObject 1
#define BCOL_set_ExtraDataInsideObject \
  BCOL_t::_3f color; \
  BCOL_t::_f transparency;
#define BCOL_set_DynamicDeltaFunction \
  ObjectData0->Velocity.y = -1;
#if set_bcol_UseEmbree
  #define BCOL_set_UseEmbree 1
#endif
#define BCOL_set_VisualSolve 1
#define BCOL_set_VisualSolve_CalculateBarycentric 1
#define BCOL_set_VisualSolve_GridContact \
  Contact = IsGridSolid(grid_result.gi);
#define BCOL_set_DPFStore \
  bcol_model_t *bcol_model; \
  uint32_t MaterialIndex; \
  fan::vec2 uv[3]; \
  fan::vec4 color;
#include <BCOL/BCOL.h>

BCOL_t* _g_bcol = new BCOL_t;
BCOL_t& g_bcol = *_g_bcol;

#include <vector>

#if set_VisualDebug == 1
  std::vector<BCOL_t::ObjectID_t> DebugBallList;
  void DebugBall(fan::vec3 pos) {
    BCOL_t::ObjectProperties_t p;
    p.Position = pos;
    BCOL_t::ObjectID_t oid = g_bcol.NewObject(&p, BCOL_t::ObjectFlag::Constant);
    DebugBallList.push_back(oid);
    g_bcol.SetObject_Velocity(oid, BCOL_t::_vf(0, 0, 0));
    g_bcol.GetObjectExtraData(oid)->color = fan::vec3(0, 1, 0);
    g_bcol.GetObjectExtraData(oid)->transparency = 0;
    BCOL_t::ShapeProperties_Circle_t sp;
    sp.Position = 0;
    sp.Size = 0.25;
    g_bcol.NewShape_Circle(oid, &sp);
  }
#endif

struct model_material_t {
  uint32_t diffuse_id;
};
std::vector<model_material_t> model_material;

#include "bcol_model.h"

enum class model_image_type : uint8_t{
  rgb,
  rgba
};

struct image_t {
  model_image_type type;
  uint8_t* data;
  fan::vec2ui size;
};

image_t model_image[999];


fan::vec4 get_color_from_image(uint32_t ti, fan::vec2 tc){
  fan::vec4 r;

  auto &image = model_image[ti];

  uint32_t cx = tc.x * image.size.x;
  uint32_t cy = tc.y * image.size.y;
  cx %= image.size.x;
  cy %= image.size.y;

  if(image.type == model_image_type::rgb){
    r[0] = f32_t(image.data[(cy * image.size.x + cx) * 3 + 0]) / 0xff;
    r[1] = f32_t(image.data[(cy * image.size.x + cx) * 3 + 1]) / 0xff;
    r[2] = f32_t(image.data[(cy * image.size.x + cx) * 3 + 2]) / 0xff;
    r[3] = 1;
  }
  else if(image.type == model_image_type::rgba){
    r[0] = f32_t(image.data[(cy * image.size.x + cx) * 4 + 0]) / 0xff;
    r[1] = f32_t(image.data[(cy * image.size.x + cx) * 4 + 1]) / 0xff;
    r[2] = f32_t(image.data[(cy * image.size.x + cx) * 4 + 2]) / 0xff;
    r[3] = f32_t(image.data[(cy * image.size.x + cx) * 4 + 3]) / 0xff;
  }
  else{
    __abort();
  }

  return r;
}

fan::vec2 CalculateBarycentric(
  const fan::vec2& b,
  const fan::vec2& v1_uv,
  const fan::vec2& v2_uv,
  const fan::vec2& v3_uv
) {
  double bz = 1.0 - b.x - b.y;
  double u = b.x * v1_uv.x + b.y * v2_uv.x + bz * v3_uv.x;
  double v = b.x * v1_uv.y + b.y * v2_uv.y + bz * v3_uv.y;
  return fan::vec2(u, v);
}


fan::vec3 reflect(const BCOL_t::_vf& incident, const BCOL_t::_vf& normal) {
  double dot_product = incident.dot(normal);
  fan::vec3 scaled_normal = normal * (2 * dot_product);
  fan::vec3 reflection = incident - scaled_normal;
  return reflection;
}

#if set_HaveGrid == 1
void PreSolve_Grid_cb(
  BCOL_t* bcol,
  const BCOL_t::ShapeInfoPack_t* sip,
  BCOL_t::_vsi32 gi,
  BCOL_t::Contact_Grid_t* c
) {
  if (IsGridSolid(gi)) {
    bcol->Contact_Grid_EnableContact(c);
    return;
  }
  bcol->Contact_Grid_DisableContact(c);
}
#endif

BCOL_t::VisualSolve_t
VisualSolve_Grid_Fragment(
  const BCOL_t::_vf& src,
  const BCOL_t::_vf& at,
  const BCOL_t::_vf& n,
  const fan::vec2& tc
) {
  BCOL_t::VisualSolve_t ret;

  __abort(); /* make stuff right */
  fan::vec3 rgb = 0;
  #if 0
  fan::vec3 rgb = get_color_from_image(0, tc);

  f32_t m = fmod(at.length(), 30) / 30;
  if (m < 0.5) {
    m = 1.0 - m;
  }
  rgb *= m * m;

  ret.at = at;
  ret.normal = reflect((at - src).normalize(), n);
  ret.multipler = 0;
  ret.rgb = rgb;
  ret.transparency = 0;
  ret.reflect = 0;
  #endif

  return ret;
}

void VisualSolve_Grid_cb(
  BCOL_t* bcol,
  BCOL_t::_vsi32 i,
  BCOL_t::_vf src,
  BCOL_t::_vf at,
  BCOL_t::VisualSolve_t* vs
) {
  BCOL_t::_vf diff = BCOL_t::_vf(i) + 0.5;
  diff -= at;

  BCOL_t::_v<BCOL_t::_dc - 1, BCOL_t::_f> barycentric = BCOL_t::CalculateBarycentricFromRectangleIntersection(BCOL_t::_vf(i) + 0.5, BCOL_t::_vf(1), at);
  BCOL_t::_vf normal = BCOL_t::CalculateNormalFromRectangleIntersection(BCOL_t::_vf(i) + 0.5, BCOL_t::_vf(1), at);
  {
    normal = BCOL_t::_vf(0);
    BCOL_t::_f ma = diff.abs().max();
    for (uint32_t d = 0; d < 3; d++) {
      if (ma == BCOL_t::abs(diff[d])) {
        normal[d] = copysign(BCOL_t::_f(1), -diff[d]);
      }
    }
  }
  *vs = VisualSolve_Grid_Fragment(src * GridBlockSize, at * GridBlockSize, normal, barycentric);
}

BCOL_t::VisualSolve_t
VisualSolve_Shape_Fragment_Circle(
  const BCOL_t::_vf& src,
  const BCOL_t::_vf& at,
  const BCOL_t::_vf& n,
  BCOL_t::_3f c,
  BCOL_t::_f transparency
) {
  BCOL_t::VisualSolve_t ret;

  ret.at = at;
  ret.normal = reflect((at - src).normalize(), n);
  ret.multipler = 0.5;
  ret.rgb = c * 0.5;
  ret.transparency = 0;
  ret.reflect = 1;

  return ret;
}

BCOL_t::VisualSolve_t
VisualSolve_Shape_Fragment_Rectangle(
  const BCOL_t::_vf& src,
  const BCOL_t::_vf& at,
  const BCOL_t::_vf& n
) {
  BCOL_t::VisualSolve_t ret;

  ret.at = at;
  ret.normal = reflect((at - src).normalize(), n);
  ret.multipler = 0.5;
  ret.rgb = fan::vec3(0.5, 0.5, 0);
  ret.transparency = 0;
  ret.reflect = 0.5;

  return ret;
}

BCOL_t::VisualSolve_t
VisualSolve_Shape_Fragment_DPF(
  const BCOL_t::_vf& src,
  const BCOL_t::_vf& at,
  const BCOL_t::_vf& n,
  BCOL_t::_v<BCOL_t::_dc - 1, BCOL_t::_f> barycentric,
  uint32_t MaterialIndex
) {
  BCOL_t::VisualSolve_t ret;

  ret.at = at;
  ret.normal = reflect((at - src).normalize(), n);
  ret.multipler = 0.75;

  fan::vec4 color = get_color_from_image(model_material[MaterialIndex].diffuse_id, barycentric);
  //fan::vec4 color = fan::vec4(fan::vec2(barycentric.x, barycentric.y), fan::vec2(1.0 - barycentric.x - barycentric.y, 1));
  ret.rgb = fan::vec3(color.x, color.y, color.z);
  //ret.rgb = fan::vec3(BCOL_t::abs(barycentric[0]), BCOL_t::abs(barycentric[1]), BCOL_t::abs(n[1]));
  //ret.rgb = fan::vec3(1, 0, 1);
  ret.transparency = 0;
  ret.reflect = 0;

  return ret;
}

void VisualSolve_Shape_cb(
  BCOL_t* bcol,
  const BCOL_t::ShapeInfoPack_t* sip,
  BCOL_t::_vf src,
  BCOL_t::_vf at,
  BCOL_t::_vf normal,
  BCOL_t::_v<BCOL_t::_dc - 1, BCOL_t::_f> barycentric,
  BCOL_t::VisualSolve_t* vs
) {
  if (sip->ShapeEnum == BCOL_t::ShapeEnum_t::Circle) {
    auto ObjectData = bcol->GetObjectData(sip->ObjectID);
    auto sd = bcol->ShapeData_Circle_Get(ObjectData->ShapeList.ptr[sip->ShapeID.ID].ShapeID);

    *vs = VisualSolve_Shape_Fragment_Circle(src, at, normal, bcol->GetObjectExtraData(sip->ObjectID)->color, bcol->GetObjectExtraData(sip->ObjectID)->transparency);
  }
  else if (sip->ShapeEnum == BCOL_t::ShapeEnum_t::Rectangle) {
    auto ObjectData = bcol->GetObjectData(sip->ObjectID);
    auto sd = bcol->ShapeData_Rectangle_Get(ObjectData->ShapeList.ptr[sip->ShapeID.ID].ShapeID);

    *vs = VisualSolve_Shape_Fragment_Rectangle(src, at, normal);
  }
  else if (sip->ShapeEnum == BCOL_t::ShapeEnum_t::DPF) {
    auto ObjectData = bcol->GetObjectData(sip->ObjectID);
    auto sd = bcol->ShapeData_DPF_Get(ObjectData->ShapeList.ptr[sip->ShapeID.ID].ShapeID);
    *vs = VisualSolve_Shape_Fragment_DPF(
      src,
      at,
      normal,
      CalculateBarycentric(
        barycentric,
        sd->u.uv[1],
        sd->u.uv[2],
        sd->u.uv[0]
      ),
      sd->u.MaterialIndex
    );
  }
  else {
    *vs = BCOL_t::VisualSolve_t(0);
  }
}

fan::vec3 rotateVector(fan::vec3 v, fan::vec2 angles) {
  f32_t y1 = v.y * std::cos(angles.x) - v.z * std::sin(angles.x);
  f32_t z1 = v.y * std::sin(angles.x) + v.z * std::cos(angles.x);

  f32_t x2 = v.x * std::cos(angles.y) + z1 * std::sin(angles.y);
  f32_t z2 = -v.x * std::sin(angles.y) + z1 * std::cos(angles.y);

  return fan::vec3(x2, y1, z2);
}

static constexpr f32_t FOV = 90.0f;

// update projection to perspective in window resize
fan::mat4 projection;
fan::mat4 view;
fan::mat4 view_projection;

fan::camera init_camera() {
  fan::camera camera;
  camera.position = { 0, 0.63, -2.96 };
  //camera.sensitivity = 0.01;
  projection = fan::mat4(1);
  view = fan::mat4(1);
  view_projection = fan::mat4(1);
  camera.update_view();
  return camera;
  //fan::vec3 camera_pos = { -2, 9, -8 };
  //fan::vec2 camera_angle = { -0.400000036, 3.29999804 };
}

fan::camera camera = init_camera();

fan::vec3 direction_vector(uint32_t xi, uint32_t yi) {

  f32_t ndcX = (2.0f * xi) / RenderSize.x - 1.0f;
  f32_t ndcY = 1.0f - (2.0f * yi) / RenderSize.y;

  fan::vec4 clipCoords = { ndcX, ndcY, -1.0f, 1.0f };
  fan::vec4 eyeCoords = view_projection.inverse() * clipCoords;

  fan::vec3 direction = fan::vec3(eyeCoords.x, eyeCoords.y, eyeCoords.z).normalize();

  return direction;
}


void ProcessSingleRay(uint8_t* Pixels, BCOL_t::_f RayPower, BCOL_t::_vf pos, BCOL_t::_vf angle, uint64_t* RayCount) {
  if (RayPower < 0.05) {
    return;
  }

  BCOL_t::VisualSolve_t vs = g_bcol.Ray(pos, angle);
  (*RayCount)++;

  vs.rgb *= RayPower;
  vs.rgb *= BCOL_t::_f(1) - (vs.reflect + vs.transparency) * vs.multipler;

  Pixels[0] += (uint32_t)(vs.rgb.x * 0xff);
  Pixels[1] += (uint32_t)(vs.rgb.y * 0xff);
  Pixels[2] += (uint32_t)(vs.rgb.z * 0xff);

  RayPower *= vs.multipler;

  ProcessSingleRay(Pixels, RayPower * vs.transparency, vs.at, angle, RayCount);
  ProcessSingleRay(Pixels, RayPower * vs.reflect, vs.at, vs.normal, RayCount);
}

uint8_t* FrameData = 0;
f64_t total_delta = 0;
BCOL_t::ObjectID_t focused_oid;

#if set_Multithread == 1
  struct RayIterateIndexes_t {
    uint32_t x0;
    uint32_t x1;
    uint32_t y0;
    uint32_t y1;
  };

  #define TRT_BME_set_Prefix fastmut
  #define TRT_BME_set_Language 1
  #define TRT_BME_set_MutexType 1
  #define TRT_BME_set_LockValue 1
  #include <WITCH/TRT/BME/BME.h>

  struct {
    RayIterateIndexes_t vec[
      (RenderSizeX / RenderSizeDivide + !!(RenderSizeX % RenderSizeDivide)) *
        (RenderSizeY / RenderSizeDivide + !!(RenderSizeY % RenderSizeDivide))
    ];
    uint32_t vec_size = 0;
    fastmut_t mut;
  }rii_pack;
#endif

void IterateRaysSingle(uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1, uint64_t* RayCount) {
  for (uint32_t y = y0; y < y1; y++) {
    for (uint32_t x = x0; x < x1; x++) {
      #if set_DisplayWindow == 0
        auto pos = BCOL_t::_vf(0, -125, 0);
        pos.x += std::sin(total_delta / 3) * 56;
        pos.z += std::cos(total_delta / 3) * 56;

        BCOL_t::_vf lookto = pos - g_bcol.GetObject_Position(focused_oid);
        lookto = lookto.square_normalize();
        fan::vec2 angle;
        angle.y = atan2(lookto.x, lookto.z);
        angle.x = -asin(lookto.y);
      #endif
      fan::vec3 angle3 = direction_vector(x, y);

      FrameData[(y * RenderSize.x + x) * 3 + 0] = 0;
      FrameData[(y * RenderSize.x + x) * 3 + 1] = 0;
      FrameData[(y * RenderSize.x + x) * 3 + 2] = 0;

      ProcessSingleRay(&FrameData[(y * RenderSize.x + x) * 3], 1, camera.position, angle3, RayCount);
    }
  }
}

struct IterateRaysMulti_RayList_t {
  uint8_t* Pixels[0x100];
  BCOL_t::_f RayPower[0x100];
  BCOL_t::_vf position[0x100];
  BCOL_t::_vf direction[0x100];
  uint32_t current = 0;
  void add(uint8_t* p_Pixels, BCOL_t::_f p_RayPower, BCOL_t::_vf p_position, BCOL_t::_vf p_direction) {
    if (p_RayPower < 0.05) {
      return;
    }
    if (current == 0x100) {
      __abort();
    }
    Pixels[current] = p_Pixels;
    RayPower[current] = p_RayPower;
    position[current] = p_position;
    direction[current] = p_direction;
    current++;
  }
  void rm(uint32_t i) {
    current--;
    Pixels[i] = Pixels[current];
    RayPower[i] = RayPower[current];
    position[i] = position[current];
    direction[i] = direction[current];
  }
};
void IterateRaysMulti_Execute(IterateRaysMulti_RayList_t& RayList) {
  BCOL_t::VisualSolve_t VisualSolve[16];
  g_bcol.Ray16(RayList.position, RayList.direction, VisualSolve);
  for (uint8_t i = 0; i < 16; i++) {
    BCOL_t::VisualSolve_t& vs = VisualSolve[i];

    vs.rgb *= RayList.RayPower[i];
    vs.rgb *= BCOL_t::_f(1) - (vs.reflect + vs.transparency) * vs.multipler;

    RayList.Pixels[i][0] += (uint32_t)(vs.rgb.x * 0xff);
    RayList.Pixels[i][1] += (uint32_t)(vs.rgb.y * 0xff);
    RayList.Pixels[i][2] += (uint32_t)(vs.rgb.z * 0xff);

    RayList.RayPower[i] *= vs.multipler;

    RayList.add(RayList.Pixels[i], RayList.RayPower[i] * vs.transparency, vs.at, RayList.direction[i]);
    RayList.add(RayList.Pixels[i], RayList.RayPower[i] * vs.reflect, vs.at, vs.normal);
  }
  for (uint8_t i = 16; i--;) {
    RayList.rm(i);
  }
}
void IterateRaysMulti(uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1, uint64_t* RayCount) {
  const uint32_t needed = 16;
  IterateRaysMulti_RayList_t RayList;

  for (uint32_t y = y0; y < y1; y++) {
    for (uint32_t x = x0; x < x1; x++) {
      BCOL_t::_vf pos = camera.position;
      #if set_DisplayWindow == 0
        pos = BCOL_t::_vf(0, -125, 0);
        pos.x += std::sin(total_delta / 3) * 56;
        pos.z += std::cos(total_delta / 3) * 56;

        BCOL_t::_vf lookto = pos - g_bcol.GetObject_Position(focused_oid);
        lookto = lookto.square_normalize();
        fan::vec2 angle;
        angle.y = atan2(lookto.x, lookto.z);
        angle.x = -asin(lookto.y);
      #endif
      fan::vec3 angle3 = direction_vector(x, y);

      FrameData[(y * RenderSize.x + x) * 3 + 0] = 0;
      FrameData[(y * RenderSize.x + x) * 3 + 1] = 0;
      FrameData[(y * RenderSize.x + x) * 3 + 2] = 0;

      RayList.add(&FrameData[(y * RenderSize.x + x) * 3], 1, pos, angle3);

      while (RayList.current >= needed) {
        IterateRaysMulti_Execute(RayList);
      }
    }
  }

  for (uint8_t i = RayList.current; i--;) {
    ProcessSingleRay(RayList.Pixels[i], RayList.RayPower[i], RayList.position[i], RayList.direction[i], RayCount);
  }
}

void IterateRays(uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1, uint64_t* RayCount) {
  #if set_bcol_IterateRaysMulti == 1
    IterateRaysMulti(x0, x1, y0, y1, RayCount);
  #else
    IterateRaysSingle(x0, x1, y0, y1, RayCount);
  #endif
}

#if set_Multithread == 1
  void GetAndIterateRays(bool* Processing, uint64_t* lock_fail_count, uint64_t* RayCount) {
    while (1) {
      RayIterateIndexes_t rii;
      #if set_Multithread_UseCond == 1
        while (rii_pack.mut.Lock()) {
          (*lock_fail_count)++;
        }
        if (rii_pack.vec_size == 0) {
          rii_pack.mut.Unlock();
          return;
        }
        *Processing = 1;
        rii = rii_pack.vec[--rii_pack.vec_size];
        rii_pack.mut.Unlock();
      #else
        while (1) {
          uint32_t size = __atomic_load_n(&rii_pack.vec_size, __ATOMIC_RELAXED);
          if (size == 0) {
            _mm_pause();
            continue;
          }

          while (rii_pack.mut.Lock()) {
            (*lock_fail_count)++;
          }

          if (rii_pack.vec_size == 0) {
            rii_pack.mut.Unlock();
            _mm_pause();
            continue;
          }

          *Processing = 1;
          rii = rii_pack.vec[--rii_pack.vec_size];
          rii_pack.mut.Unlock();
          break;
        }
      #endif
      IterateRays(rii.x0, rii.x1, rii.y0, rii.y1, RayCount);
      *Processing = 0;
    }
  }

  void GetAndIterateRays_Main(uint64_t* lock_fail_count, uint64_t* RayCount) {
    #if set_Multithread_UseCond == 1
    bool Processing = 0;
    GetAndIterateRays(&Processing, lock_fail_count, RayCount);
    #else
    while (1) {
      while (rii_pack.mut.Lock()) {
        (*lock_fail_count)++;
      }
      if (rii_pack.vec_size == 0) {
        rii_pack.mut.Unlock();
        return;
      }
      auto rii = rii_pack.vec[--rii_pack.vec_size];
      rii_pack.mut.Unlock();
      IterateRays(rii.x0, rii.x1, rii.y0, rii.y1, RayCount);
    }
    #endif
  }

  struct tp_pack_t {
    #if set_Multithread_UseCond == 1
      #define TRT_BME_set_Prefix Cond
      #define TRT_BME_set_Language 1
      #define TRT_BME_set_AreWeInsideStruct 1
      #define TRT_BME_set_Conditional
      #include <WITCH/TRT/BME/BME.h>
      Cond_t cond;
      bool ping = 0;
    #endif
    bool Processing = 0;
    uint64_t rii_lock_fail_count = 0;
    uint64_t RayCount = 0;
  };

  void cb_thread(void* p) {
    tp_pack_t* tp_pack = (tp_pack_t*)p;
    while (1) {
      #if set_Multithread_UseCond == 1
        tp_pack->cond.Lock();
        while (tp_pack->ping == 0) {
          tp_pack->cond.Wait();
        }
        tp_pack->ping = 0;
        tp_pack->cond.Unlock();
      #endif

      GetAndIterateRays(&tp_pack->Processing, &tp_pack->rii_lock_fail_count, &tp_pack->RayCount);
    }
  }

  uint64_t Main_rii_lock_fail_count = 0;
  uint64_t sleep_wait_count = 0;
#endif

uint64_t MainRayCount = 0;

int main() {
  #if set_Multithread == 1
    uint32_t thread_count = WITCH_num_online_cpus();
    if (thread_count > 0x100) {
      fprint("(thread_count > 0x100)");
      return 0;
    }
    tp_pack_t* tp_pack = new tp_pack_t[thread_count - 1];

    for (uint32_t i = 0; i < (thread_count - 1); i++) {
      TH_open((void*)cb_thread, &tp_pack[i]);
    }
  #endif

  {
    BCOL_t::OpenProperties_t op;
    #if set_HaveGrid == 1
      op.GridBlockSize = GridBlockSize;
    #endif
    g_bcol.Open(op);
    #if set_HaveGrid == 1
      g_bcol.PreSolve_Grid_cb = PreSolve_Grid_cb;
      g_bcol.VisualSolve_Grid_cb = VisualSolve_Grid_cb;
    #endif
    g_bcol.VisualSolve_Shape_cb = VisualSolve_Shape_cb;
    g_bcol.PreSolve_Shape_cb = [](
      BCOL_t* bcol,
      const BCOL_t::ShapeInfoPack_t* sip0,
      const BCOL_t::ShapeInfoPack_t* sip1,
      BCOL_t::Contact_Shape_t* c
      ) {
        bcol->Contact_Shape_EnableContact(c);
      };
  }

  bcol_model_t bcol_model("player.gltf");
  bcol_model.open();
  auto found = bcol_model.animation_list.find("Idle");
  if (found != bcol_model.animation_list.end()) {
    found->second.weight = 1.0;
  }

  {
    uint32_t i = 0;
    for(auto & t : fan_3d::model::cached_texture_data){
      model_image[i].data = t.second.data.data();
      model_image[i].size = t.second.size;
      if(t.second.channels == 3){
        model_image[i].type = model_image_type::rgb;
      }
      else if(t.second.channels == 4){
        model_image[i].type = model_image_type::rgba;
      }
      else{
        __abort();
      }
      i++;
    }
  }

  //bcol_model.fms.m_transform = bcol_model.fms.m_transform.scale(1);

  //bcol_model.fms.m_transform = bcol_model.fms.m_transform.rotate(fan::math::radians(180.f), fan::vec3(0, 1, 0));
  //
  {
    // auto anid = fms.create_an("an_name", weight, duration);
    //auto animation_node_id = fms.fk_set_rot(anid, "Left_leg", 0.001/* time in seconds */,
    //  fan::vec3(1, 0, 0), 0
    //);
    //auto animation_node_id3 = fms.fk_set_rot(anid, "Left_leg", 1/* time in seconds */,
    //  fan::vec3(1, 1, 0), fan::math::pi
    //);

    //auto anid = bcol_model.fms.create_an("an_name", 0.5);
    ////auto anid2 = bcol_model.fms.create_an("an_name2", 0.5);

    //auto animation_node_id1 = bcol_model.fms.fk_set_rot(anid, "Armature_Chest", 0.001/* time in seconds */,
    //  fan::vec3(1, 0, 0), 0
    //);
    //auto animation_node_id = bcol_model.fms.fk_set_rot(anid, "Armature_Chest", 0.3/* time in seconds */,
    //  fan::vec3(1, 0, 0), fan::math::pi / 2
    //);

    //auto animation_node_id3 = bcol_model.fms.fk_set_rot(anid, "Armature_Chest", 0.6/* time in seconds */,
    //  fan::vec3(1, 0, 0), -fan::math::pi / 3
    //);

    //auto animation_node_id2 = bcol_model.fms.fk_set_rot(anid, "Armature_Upper_Leg_L", 0.6/* time in seconds */,
    //  fan::vec3(1, 0, 0), -fan::math::pi
    //);
    //auto animation_node_id2 = bcol_model.fms.fk_set_rot(anid, 0.4/* time in seconds */, "Armature_Lower_Arm_R", fan::vec3(0, 180, 0));

    //auto animation_node_id3 = bcol_model.fms.fk_set_rot(anid2, 0.5/* time in seconds */, "Armature_Lower_Leg_L", fan::vec3(0, 180, 0));
    //auto animation_node_id4 = bcol_model.fms.fk_set_rot(anid2, 0.7/* time in seconds */, "Armature_Lower_Leg_L", fan::vec3(0, 180, 0));
  }

  uint64_t FrameCount = 0;

  FrameData = (uint8_t*)malloc(RenderSize.x * RenderSize.y * 3);

  #if set_DisplayWindow == 1
    loco_t loco = loco_t::properties_t{ .window_size = {RenderSize.x, RenderSize.y} };
    loco.camera_set_ortho(loco.orthographic_camera.camera, 
      fan::vec2(-1, 1),
      fan::vec2(-1, 1)
    );

    projection = fan::math::perspective<fan::mat4>(fan::math::radians(FOV), (f32_t)gloco->window.get_size().x / (f32_t)gloco->window.get_size().y, 0.1f, 1000.0f);

    loco.window.add_mouse_motion([&](const auto& d) {
      if (ImGui::IsMouseDown(ImGuiMouseButton_Middle)) {
        camera.rotate_camera(d.motion);
      }
    });

    loco_t::image_t image;
    image = loco.image_create();
    loco_t::image_load_properties_t lp;
    lp.format = fan::opengl::GL_RGB;
    lp.internal_format = fan::opengl::GL_RGB;

    loco_t::sprite_t::properties_t p;
    p.size = fan::vec2(1, 1);
    p.position = fan::vec3(0, 0, 0x7f);
    p.image = image;
    loco_t::shape_t id = p;

    loco.set_vsync(false);
    loco.console.commands.call("show_fps 1");
  #endif

  uint64_t RayCount = 0;
  uint64_t frame = 0;

  #if set_DisplayWindow == 1
    loco.loop([&]
  #else
    while (1)
  #endif
  {
      if (frame == 1) {
        loco.console.commands.call("set_target_fps 0");
        
      }
      frame++;
    #if set_DisplayWindow == 0
      if (FrameCount >= set_RenderFrameCount) {
        return 0;
      }
    #endif

    #if set_DisplayWindow == 1
      if (!image.iic()) {
        loco.image_unload(image);
      }
    #endif

    #if set_VisualDebug == 1
      for (uintptr_t i = 0; i < DebugBallList.size(); i++) {
        g_bcol.UnlinkObject(DebugBallList[i]);
        g_bcol.RecycleObject(DebugBallList[i]);
      }
      DebugBallList.clear();
    #endif

    bcol_model.delta(0.01);
    g_bcol.Step(0.01);
    g_bcol.BakeCurrentForVisualSolve();

    f64_t t0 = T_nowf();
    #if set_Multithread == 1
      #if set_Multithread_UseCond == 1
        while (rii_pack.mut.Lock()) {
          Main_rii_lock_fail_count++;
        }
      #endif
      for (uint32_t y = 0; y < RenderSize.y;) {
        uint32_t lefty = RenderSizeY - y;
        if (lefty > RenderSizeDivide) {
          lefty = RenderSizeDivide;
        }
        #if set_Multithread_UseCond == 0
        while (rii_pack.mut.Lock()) {
          Main_rii_lock_fail_count++;
        }
        #endif
        for (uint32_t x = 0; x < RenderSize.x;) {
          uint32_t leftx = RenderSizeX - x;
          if (leftx > RenderSizeDivide) {
            leftx = RenderSizeDivide;
          }
          rii_pack.vec[rii_pack.vec_size++] = { .x0 = x, .x1 = x + leftx, .y0 = y, .y1 = y + lefty };
          x += leftx;
        }
        #if set_Multithread_UseCond == 0
        rii_pack.mut.Unlock();
        #endif
        y += lefty;
      }

      #if set_Multithread_UseCond == 1
        rii_pack.mut.Unlock();
      #endif

      #if set_Multithread_UseCond == 1
        for (uint32_t ti = 0; ti < thread_count - 1; ti++) {
          tp_pack[ti].cond.Lock();
          tp_pack[ti].ping = 1;
          tp_pack[ti].cond.Signal();
          tp_pack[ti].cond.Unlock();
        }
      #endif

      GetAndIterateRays_Main(&Main_rii_lock_fail_count, &MainRayCount);
    #else
      IterateRays(0, RenderSize.x, 0, RenderSize.y, &MainRayCount);
    #endif

    #if set_PrintStats == 1
      #if set_Multithread == 0
        fprint("frame", FrameCount);
      #else
        fprint("frame", FrameCount, "locks", Main_rii_lock_fail_count, sleep_wait_count);
        for (uint32_t i = 0; i < thread_count - 1; i++) {
          uint64_t rii_lock_fail_count = __atomic_load_n(&tp_pack[i].rii_lock_fail_count, __ATOMIC_RELAXED);
          fprint("thread", i + 1, rii_lock_fail_count);
        }
      #endif
      fprint("total rays >=", RayCount);
    #endif

    RayCount = MainRayCount;

    #if set_Multithread == 1
      for (uint32_t i = 0; i < thread_count - 1;) {
        RayCount += __atomic_load_n(&tp_pack[i].RayCount, __ATOMIC_RELAXED);

        while (1) {
          bool Processing = __atomic_load_n(&tp_pack[i].Processing, __ATOMIC_RELAXED);
          if(Processing == 0){
            i++;
            break;
          }
          sleep_wait_count++;
          _mm_pause();
        }
      }
    #endif

    f64_t t1 = T_nowf();
    //fprint("t is", t1 - t0);

    #if set_stdout_data == 1
      print("%.*s", (uintptr_t)RenderSize.x * RenderSize.y * 3, FrameData);
    #endif

    #if set_DisplayWindow == 1
      fan::image::image_info_t ii;
      ii.data = (uint8_t*)FrameData;
      ii.size.x = RenderSize.x;
      ii.size.y = RenderSize.y;
      if (image.iic() == false) {
        loco.image_unload(image);
        image.sic();
      }
      image = loco.image_load(ii, lp);
    #endif
    #if set_DisplayWindow == 1
      if(gloco->window.key_pressed(fan::key_left_alt)){
        camera.move(1);
      }
      else {
        camera.move(100);
      }

      view = fan::math::look_at_left<fan::mat4, fan::vec3>(0, camera.m_front, camera.m_up);
      view_projection = projection * view;

      //fan::ray3_t ray = gloco->convert_mouse_to_ray(camera.position, projection, view);
      //bcol_model.fms.iterate_joints(bcol_model.fms.parsed_model.model_data.skeleton, [&](fan_3d::animation::joint_t& joint){
      //  fan::vec3 position = joint.global_transform.get_translation();
      //  fan::vec3 size = 1;
      //  if (gloco->is_ray_intersecting_cube(ray, position, size)) {
      //    //fan::print("a");
      //  }
      //});

      total_delta += loco.delta_time;
    #else
      total_delta += set_RenderTick;
    #endif

    FrameCount++;
  }
  #if set_DisplayWindow == 1
    );
  #endif

  return 0;
}
