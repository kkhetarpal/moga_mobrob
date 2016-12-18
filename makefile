default: 
	g++ -std=c++0x -ggdb -DGALLDATA $(ap_flags) \
`pkg-config --cflags playerc++` \
`pkg-config --libs playerc++`  \
-I/usr/local/include/player-2.1 \
-L/usr/lib  -lgsl -lgslcblas \
*.c fuzzy.cpp simulator.cpp ap_nsga2_utils.cpp mlp.cpp

drive: 
	g++ -std=c++0x -ggdb -DGALLDATA $(ap_flags) \
`pkg-config --cflags playerc++` \
`pkg-config --libs playerc++`  \
-I/usr/local/include/player-2.1 \
-L/usr/lib  -lgsl -lgslcblas \
drive.cpp

test_simulator: 
	g++ -std=c++0x -ggdb $(ap_flags) \
`pkg-config --cflags playerc++` \
`pkg-config --libs playerc++`  \
-I/usr/local/include/player-2.1 \
-L/usr/lib  -lgsl -lgslcblas \
test_simulator.cpp simulator.cpp fuzzy.cpp mlp.cpp

tags:
	ctags -e -R --c++-kinds=+p --c-kinds=+p

clean:
	rm -f a.out *.o *.txt *.log
