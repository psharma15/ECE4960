import numpy as np
from decimal import Decimal, getcontext
import logging
import warnings
from math import factorial

getcontext().prec = 40
	
def integerOverflow():
	intFactorial = np.int32(1)
	for i in range(1,15):
		try:
			intFactorial = ((intFactorial) * i)
			print(i, intFactorial)
		except RuntimeWarning as e:
			  file = open('EngineerReport.txt','w')
			  file.write("############################\n# INTEGER OVERFLOW TESTING #\n############################\n")
			  file.write("Runtime Warning Error below occured at iteration %d\n##%s##\n\n"%(i,e))
			  file.close()
			  print(i, intFactorial)	
	print ("\n")
	

def integerDivideZero():
	num1 = 5
	zero = 0
	try:
		val = num1/zero
	except ZeroDivisionError as  e:
		file = open('EngineerReport.txt','a')
		file.write("####################################\n# INTEGER DIVISION BY ZERO TESTING #\n####################################\n")
		file.write("The ZeroDivisionError below occured \n##%s##\n\n"%e)
		file.close()
	print ("\n")
	

def floatOverflow():
	float1 = np.float64(1.02)

	for i in range(1,200):
		try:
			float1 = float1 * i
			print ("%d, %.30f" %(i, float1))
		except RuntimeWarning as e:
			file = open('EngineerReport.txt','a')
			file.write("###########################\n# FLOAT OVERFLOW TESTING #\n###########################\n")
			file.write("Runtime Warning Error below occured at iteration %d\n##%s##\n\n"%(i,e))
			file.close()
			print ("%d, %.30f" %(i, float1))
			break
	print ("\n")
	

def floatOpINF():
	val1 = np.inf
	val2 = -np.inf
	
	try:
		test1 = np.sin(val1)
	except RuntimeWarning as e:
		file = open('EngineerReport.txt','a')
		file.write("###############################\n# FLOAT INF AND NINF TESTING #\n###############################\n")
		file.write("Runtime Warning Error below occured \n##%s##\n\n"%(e))
		
	test2 = 1/ val1
	test3 = np.exp(val1)

	try:
		test4 = np.sin(val2)

	except RuntimeWarning as e:
		file.write("Runtime Warning Error below occured \n##%s##\n\n"%(e))
		file.close()	
	
	test5 = 1/val2
	test6 = np.exp(val2)
	
	print (test1)
	print (test2)
	print (test3)
	print (test4)
	print (test5)
	print (test6)


	#Testing the propagation
	print(np.inf + np.inf)
	print(0*np.inf)
	print(np.inf - np.inf)
	print ("\n")


def floatOpNAN():
	#nan1 = np.nan
	#nan1 = np.float64(np.sqrt(-1.0))
	nan1 = (np.inf/np.inf)
	print (nan1)
	print (np.isnan(nan1))

	#Testing the propagation
	print(1/nan1)
	print(nan1 + np.inf)
	print(0 * nan1)
	print ("\n")



def signedZero():
	pZero = 1.0/(np.inf)
	nZero = 1.0/(-np.inf)

	psmallNo = 1e-323
	nsmallNo = -1e-323
	file = open('EngineerReport.txt','a')
	file.write("###############################\n# SIGNED ZERO TESTING #\n###############################\n")
	try:
		plog = np.log(pZero)
	except RuntimeWarning as e:
		file.write("Runtime Warning Error below occured \n##%s##\n\n"%(e))

	try:
		nlog = np.log(nZero)
	except RuntimeWarning as e:
		file.write("Runtime Warning Error below occured \n##%s##\n\n"%(e))	
	
	pSin = (np.sin(psmallNo))/psmallNo
	nSin = (np.sin(nsmallNo))/nsmallNo

	aSin = (np.sin(nsmallNo))/abs(nsmallNo)

	print (plog)
	print (nlog)
	print (pSin)
	print (nSin)
	print (aSin)
	print ("\n")


def floatSoftLand():
	a = 308
	b = 315
	c = 1.2345678912345*(10**(-a))
	
	for  i in range(20):
		x = (10**(-a))
		print("sin part is %.14f"%((np.sin(1.23456789012345*x))/x))
		a = a + 1
		y = 1.23456789 *(10**(-b))
		print("\nx - y = %.14f\n" %(c-y))
		
		try:
			print("\n x/y = %.14f" %(c/y))
		except ZeroDivisionError as e:
			pass

		b = b + 1



def piCalculator():
	prec = 25
	a = Decimal(0) 
	
	for k in range(prec):
		b = ((-1)**k)*(factorial(6*k)) * (13591409 +(545140134*k))
		c = (factorial(3*k))*((factorial(k))**3)
		d = 640320 **((3*k) + Decimal(1.5))
		a = Decimal(a) +((Decimal(b))/(Decimal((Decimal(c))*(Decimal(d)))))
	final = Decimal(12) * Decimal(a)	
	pi = Decimal(1)/Decimal(final)

	return pi


integerOverflow()
integerDivideZero()
floatOverflow()
floatOpINF()
floatOpNAN()
signedZero()
floatSoftLand()
print (piCalculator())