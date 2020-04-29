CC := g++
CFLAGS := -Wall -g -Isrc
TARGET := rosenbrock_2d

.PHONY: tests clean plot

tests: $(TARGET)

$(TARGET):
	@$(CC) $(CFLAGS) test/$@.cpp -o bin/$@
	@bin/$@

clean:
	rm -rf bin/*

plot:
	gnuplot -e "load 'img/viridis.pal'; splot '$(file)' with pm3d; pause -1"