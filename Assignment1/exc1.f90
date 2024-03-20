program exc1
	implicit none
		!-32768 to 32767
	INTEGER*4 :: a, b, res1
	double precision :: PI, pow1, pow2, sqrr1, res2, res3, res4
	pow1 = 1d32
	pow2 = 1d21
	sqrr1 = sqrt(2d0)
	PI = 4.D0*DATAN(1.D0)
	a = 2000000
	b = 1

	res1 = a+b

	res2 = PI*pow1
	res3 = sqrr1*pow2
	res4 = res2+res3

	print *,'EXC A:', a, b, 'the sum is', res1
	print *,'EXC B_1:', PI, pow1, 'the multiplication is', res2
	print *,'EXC B_2:', sqrr1, pow2, 'the multiplication is', res3
	print *, '----------------'
	print *, 'EXC B final result is', res4
end program exc1
