test.exe: testsrc/test.o cpp/gaussian_rand.o cpp/Particle.o cpp/System.o cpp/MT.o
	icc -o $@ $^

BD.exe: src/main.o cpp/gaussian_rand.o cpp/Particle.o cpp/System.o cpp/MT.o
	icc -o $@ $^

measuretime.exe: measuresrc/timer.o cpp/gaussian_rand.o cpp/Particle.o cpp/System.o cpp/MT.o
	icc -o $@ $^

%.o: %.cpp h/%.hpp h/parameters.hpp
	icc -o $@ -c $<

%.o: %.cpp h/parameters.cuh
	icc -o $@ -c $<