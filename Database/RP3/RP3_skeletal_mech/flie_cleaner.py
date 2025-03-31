import os

file_ = "mechanism.inp"

with open(file_,"r") as read:
	mech = read.readlines()
	
string = ""
for line in mech:
	if "!" in line:
		continue
	string+=line

clean_mech = open("clean_mech.inp","w").write(string)
