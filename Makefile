CC = gcc
CFLAGS = -c -pedantic -Wall
LIBS = -lm # Dodanie bibliotek statycznych i dzielonych (w tym przypadku potrzebnajest np. biblioteka math.h)
STANDARD = -std=c99

all: dsplib

dsplib: main.o dsplib.o
	@echo Konsolidacja 
	$(CC) main.o dsplib.o $(STANDARD) $(LIBS) -o dsplib

main.o: main.c
	@echo Tworzenie obiektu main.o
	$(CC) $(CFLAGS) main.c

dsplib.o: dsplib.c dsplib.h
	@echo Tworzenie obiektu dsplib.o
	$(CC) $(CFLAGS)  dsplib.c

clean:
	@echo Usuwanie wszystkich obiektow
	rm -rf main.o dsplib.o
