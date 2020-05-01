CC := g++
CFLAGS := -Wall -g -Isrc -DDEBUG
TARGET := ackley2d rosenbrock_2d

.PHONY: test clean plot $(TARGET)

build: $(TARGET)

$(TARGET):
	@$(CC) $(CFLAGS) test/$@.cpp -o bin/$@

clean:
	rm -rf bin/*

plot:
ifeq ($(strip $(fun)),)
	@printf "function file not set !\n"
else ifeq ($(strip $(path)),)
	gnuplot -e "load 'img/viridis.pal'; set contour base; set cntrlabel onecolor; set cntrparam levels 10; splot '$(fun)' with pm3d; pause -1"
else
	gnuplot -e "set style line 10 lt rgb 'white' lw 2 pt 2; load 'img/viridis.pal'; set contour base; set cntrlabel onecolor; set cntrparam levels 10; splot '$(fun)' with pm3d, '$(path)'  with linespoints ls 10; pause -1"
endif