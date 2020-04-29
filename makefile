CC := g++
CFLAGS := -Wall -g -Isrc
TARGET := rosenbrock_2d

.PHONY: tests clean

tests: $(TARGET)

$(TARGET):
	@$(CC) $(CFLAGS) test/$@.cpp -o bin/$@
	@bin/$@

clean:
	rm -rf bin/*