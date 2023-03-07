#!/bin/bash

make clean;
make;
./main;
python vis.py;
ffmpeg -y -i anim/%d.png -framerate 480 out1.mp4;
rm anim/* 
