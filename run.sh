
g++ -o burgers burgers.cc
./burgers
cd ./output
 convert -layers optimize -loop 0 -delay 50 0000.jpg -delay 10 *.jpg anime.gif
cd ../
eog ./output/anime.gif
