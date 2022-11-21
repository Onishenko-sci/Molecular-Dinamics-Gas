all:
	g++ main2.cpp ./code/Vector2.cpp ./code/molecular_dinamics.cpp -o MakeMD.out
	./MakeMD.out
	python3 ShowMD.py