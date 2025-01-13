#ifndef RenderSizeX
  #define RenderSizeX 512
#endif
#ifndef RenderSizeY
  #define RenderSizeY 512
#endif

#include <WITCH/WITCH.h>
#include <WITCH/IO/IO.h>
#include <WITCH/IO/print.h>
#include <WITCH/A/A.h>

#include <termios.h>

void _print(int fd, const char* format, ...){
  IO_fd_t iofd;
  IO_fd_set(&iofd, fd);
  va_list argv;
  va_start(argv, format);
  IO_vprint(&iofd, format, argv);
  va_end(argv);
}
#define printout(...) _print(STDOUT_FILENO, __VA_ARGS__)
#define printerr(...) _print(STDERR_FILENO, __VA_ARGS__)

constexpr f32_t GridBlockSize = 64;

#include <fan/types/types.h>
#include <fan/types/vector.h>

#define BCOL_set_prefix BCOL
#define BCOL_set_SupportGrid 1
#define BCOL_set_Dimension 2
#define BCOL_set_DynamicDeltaFunction \
  ObjectData0->Velocity[0] += 0; \
  ObjectData0->Velocity[1] += 200 * delta;
#include <BCOL/BCOL.h>

bool IsGridSolid(sint32_t *p){
  if(p[0] == 0 && p[1] == 2){
    return true;
  }
  for(uint32_t grid_d = 0; grid_d < BCOL_t::_dc; grid_d++){
    if(p[grid_d] > 2 || p[grid_d] < -3){
      return true;
    }
  }
  return false;
}

void PreSolve_Grid_cb(
  BCOL_t* bcol,
  const BCOL_t::ShapeInfoPack_t* sip,
  BCOL_t::_vsi32 gi,
  BCOL_t::Contact_Grid_t* c
){
  if(IsGridSolid((sint32_t *)&gi)){
    bcol->Contact_Grid_EnableContact(c);
    return;
  }
  bcol->Contact_Grid_DisableContact(c);
}

int main(){
  {
    int flags = fcntl(STDIN_FILENO, F_GETFL, 0);
    if (flags == -1) {
      perror("fcntl(F_GETFL)");
      return -1;
    }

    flags |= O_NONBLOCK;

    if (fcntl(STDIN_FILENO, F_SETFL, flags) == -1) {
      perror("fcntl(F_SETFL)");
      return -1;
    }
  }
  {
    struct termios term;

    tcgetattr(STDIN_FILENO, &term);

    term.c_lflag &= ~(ICANON | ECHO | ISIG);
    term.c_iflag &= ~(IXON | ICRNL | INLCR);
    term.c_oflag &= ~OPOST;

    tcsetattr(STDIN_FILENO, TCSANOW, &term);
  }

  BCOL_t g_bcol;
  {
    BCOL_t::OpenProperties_t op;
    op.GridBlockSize = GridBlockSize;
    g_bcol.Open(op);
  }
  g_bcol.PreSolve_Grid_cb = PreSolve_Grid_cb;
  g_bcol.PreSolve_Shape_cb = [](
    BCOL_t* bcol,
    const BCOL_t::ShapeInfoPack_t* sip0,
    const BCOL_t::ShapeInfoPack_t* sip1,
    BCOL_t::Contact_Shape_t* c
    ){
      bcol->Contact_Shape_EnableContact(c);
    };

  uint8_t *FrameData = (uint8_t*)malloc(RenderSizeX * RenderSizeY * 3);

  BCOL_t::ObjectProperties_t op0;
  op0.Position = 0;
  op0.Position[0] = -64;
  op0.Position[1] = +64;
  auto SObjectID = g_bcol.NewObject(&op0, BCOL_t::ObjectFlag::Constant);
  BCOL_t::ShapeProperties_Rectangle_t sp0;
  sp0.Position = 0;
  sp0.Size = 48;
  sp0.Size[0] = 48;
  //sp0.Rotation = BCOL_t::_rotf(2);
  sp0.Rotation = 0;
  auto SShapeID = g_bcol.NewShape_Rectangle(SObjectID, &sp0);

  BCOL_t::ObjectProperties_t op;
  op.Position = 0;
  op.Position[0] = -112 - 4;
  op.Position[1] = 112 + 48;
  auto ObjectID = g_bcol.NewObject(&op, 0);
  BCOL_t::ShapeProperties_Circle_t sp;
  sp.Position = 0;
  sp.Size = 2;
  auto ShapeID = g_bcol.NewShape_Circle(ObjectID, &sp);

  BCOL_t::_vf campos = 0;
  BCOL_t::_vf fpixelpos = 0;

  while(1){
    {
      IO_fd_t iofd;
      IO_fd_set(&iofd, STDIN_FILENO);
      uint8_t data[0x1000];
      auto read_size = IO_read(&iofd, data, sizeof(data));
      if(read_size < 0){
        __abort();
      }
      for(uintptr_t i = 0; i < read_size; i++){
        BCOL_t::_vf vel = g_bcol.GetObject_Velocity(ObjectID);
        if(data[i] == 'a'){
          vel[0] -= 50;
        }
        else if(data[i] == 's'){
          vel[1] += 50;
        }
        else if(data[i] == 'd'){
          vel[0] += 50;
        }
        else if(data[i] == 'w'){
          vel[1] -= 250;
        }
        else if(data[i] == 'q'){
          if constexpr(BCOL_t::_dc >= 3){
            fpixelpos[2] += 1;
          }
        }
        else if(data[i] == 'e'){
          if constexpr(BCOL_t::_dc >= 3){
            fpixelpos[2] -= 1;
          }
        }
        else{
          goto gt_end_loop;
        }
        g_bcol.SetObject_Velocity(ObjectID, vel);
      }
    }
    g_bcol.Step(0.01);

    for(uint32_t y = 0; y < RenderSizeY; y++){
      for(uint32_t x = 0; x < RenderSizeX; x++){
        fan::vec2 _fpixelpos{x, y};
        _fpixelpos -= fan::vec2(RenderSizeX, RenderSizeY) / 2;
        fpixelpos[0] = _fpixelpos[0];
        fpixelpos[1] = _fpixelpos[1];
        BCOL_t::_vsi32 gridpos = fpixelpos / GridBlockSize;
        for(auto dc = BCOL_t::_dc; dc--;){
          if(fpixelpos[dc] < 0){
            gridpos[dc]--;
          }
        }
        bool igs = IsGridSolid((sint32_t *)&gridpos);
        FrameData[(y * RenderSizeX + x) * 3 + 0] = (uint8_t)igs * 0xff;
        FrameData[(y * RenderSizeX + x) * 3 + 1] = (uint8_t)igs * 0xff;
        FrameData[(y * RenderSizeX + x) * 3 + 2] = (uint8_t)igs * 0xff;

        auto pObjectID = g_bcol.GetObjectIDByPosition(fpixelpos);
        if(pObjectID == ObjectID){
          FrameData[(y * RenderSizeX + x) * 3 + 0] = 0xff;
          FrameData[(y * RenderSizeX + x) * 3 + 1] = 0;
          FrameData[(y * RenderSizeX + x) * 3 + 2] = 0;
        }
        else if(!pObjectID.iic()){
          FrameData[(y * RenderSizeX + x) * 3 + 0] = 0;
          FrameData[(y * RenderSizeX + x) * 3 + 1] = 0xff;
          FrameData[(y * RenderSizeX + x) * 3 + 2] = 0;
        }
      }
    }

    printout("%.*s", (uintptr_t)RenderSizeX * RenderSizeY * 3, FrameData);

    auto sd = g_bcol.ShapeData_Rectangle_Get(SShapeID);
    sd->Rotation += fan::math::pi;
  }

  gt_end_loop:;

  {
    struct termios term;

    tcgetattr(STDIN_FILENO, &term);
    term.c_lflag |= (ICANON | ECHO | ISIG);
    term.c_iflag |= (IXON | ICRNL | INLCR);
    term.c_oflag |= OPOST;
    tcsetattr(STDIN_FILENO, TCSANOW, &term);
  }

  {
    IO_fd_t iofd;
    IO_fd_set(&iofd, STDOUT_FILENO);
    IO_close(&iofd);
  }

  return 0;
}
