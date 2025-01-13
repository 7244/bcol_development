clang++ -std=c++20 main.cpp -o a.exe -s -O3 -march=native -mtune=native

./a.exe | ffplay -f rawvideo -loglevel quiet -video_size 512x512 -framerate 60 -pixel_format rgb24 -
