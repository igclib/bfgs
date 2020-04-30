CC := g++
CFLAGS := -Wall -g -Isrc -DDEBUG
TARGET := rosenbrock_2d

.PHONY: test clean plot

build: $(TARGET)

test:$(TARGET)
	bin/$^

$(TARGET):
	$(CC) $(CFLAGS) test/$@.cpp -o bin/$@
	

clean:
	rm -rf bin/*

plot:
ifeq ($(strip $(file)),)
	@printf "file not set !\n"
else
	gnuplot -e "load 'img/viridis.pal'; splot '$(file)' with pm3d; pause -1"
endif