all:
	g++ main.cpp ./code/Vector2.cpp -o MakeMD.out
	./MakeMD.out
	python3 ShowMD.py