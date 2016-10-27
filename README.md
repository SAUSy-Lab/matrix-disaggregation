# matrix-disaggregation
synthetically disaggregate flows for flow matrix with disaggregated O's &amp; D's plus OD costs

# Written by WRC on 05/12/2014

a.To Compile the code
  make

b.There is a test example (Wichita) in `pwd`/example/1
  In each of the the original input files the first line was removed 
  so that I can be used by my C code.
  To do the test please ./cscript  in the directory `pwd`/example/1

  If you want to see all the input options, just type ./geo

I also wrote some simple Python code which allows you to compare the 
original basic code output with the output of the new C-code

The python script can be found in  `pwd`/example/1 and it is 
advised to use python 2.7x or python 3.x
When you are in `pwd`/example/1, type
python check.py 
