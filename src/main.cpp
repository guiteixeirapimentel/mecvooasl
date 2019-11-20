#include <iostream>
#include <vector>
#define PI 3.14159265359

#include "f1.h"
#include "f2.h"
#include "f3.h"
#include "f4.h"

#include "Matriz.h"
#include "MathUtils.h"

int main()
{
	const double S = 16.0;
	const double meanChord = 1.5;
	const double mass = 1250.0;
	const double momentInertia = 1825.0;
	const double rho = 1.0;
	const double gravity = 9.80665;

	DragCoefficient CD(0.02806, 0.121);
	LiftCoefficient CL(0.2977, 4.41, 1.7, 3.9, 0.43);
	MomentCoefficient CM(0.03917, -0.613, -7.27, -12.4, -1.122);

	XForce xf(CL, CD);
	ZForce zf(CL, CD);
	YMoment ym(CM);

	ThrustForce T(1100.0);

	f1 errorFunc1(xf, T, rho, gravity, S);
	f2 errorFunc2(zf, rho, gravity, S);
	f3 errorFunc3(ym, rho, gravity, S);
	f4 errorFunc4;

	State Xref(70.0, 0.0, 0.0, 0.0);
	StateDot Xdotref(0.0, 0.0, 0.0, 0.0);
	double deltaEref = 2 * PI / 180.0;

	std::vector<double> EE;
	EE.resize(4*4);

	EE[0] = errorFunc1.JacobianX1(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[1] = errorFunc1.JacobianX2(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[2] = errorFunc1.JacobianX3(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[3] = errorFunc1.JacobianX4(meanChord, Xref, Xdotref, mass, deltaEref);

	EE[0 + (4 * 1)] = errorFunc2.JacobianX1(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[1 + (4 * 1)] = errorFunc2.JacobianX2(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[2 + (4 * 1)] = errorFunc2.JacobianX3(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[3 + (4 * 1)] = errorFunc2.JacobianX4(meanChord, Xref, Xdotref, mass, deltaEref);

	EE[0 + (4 * 2)] = errorFunc3.JacobianX1(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[1 + (4 * 2)] = errorFunc3.JacobianX2(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[2 + (4 * 2)] = errorFunc3.JacobianX3(meanChord, Xref, Xdotref, mass, deltaEref);
	EE[3 + (4 * 2)] = errorFunc3.JacobianX4(meanChord, Xref, Xdotref, mass, deltaEref);

	EE[0 + (4 * 3)] = errorFunc4.JacobianX1(Xref, Xdotref);
	EE[1 + (4 * 3)] = errorFunc4.JacobianX2(Xref, Xdotref);
	EE[2 + (4 * 3)] = errorFunc4.JacobianX3(Xref, Xdotref);
	EE[3 + (4 * 3)] = errorFunc4.JacobianX4(Xref, Xdotref);

	std::vector<double> AA;
	AA.resize(4*4);

	AA[0] = errorFunc1.JacobianXdot1(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[1] = errorFunc1.JacobianXdot2(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[2] = errorFunc1.JacobianXdot3(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[3] = errorFunc1.JacobianXdot4(meanChord, Xref, Xdotref, mass, deltaEref);

	AA[0 + (4 * 1)] = errorFunc2.JacobianXdot1(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[1 + (4 * 1)] = errorFunc2.JacobianXdot2(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[2 + (4 * 1)] = errorFunc2.JacobianXdot3(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[3 + (4 * 1)] = errorFunc2.JacobianXdot4(meanChord, Xref, Xdotref, mass, deltaEref);

	AA[0 + (4 * 2)] = errorFunc3.JacobianXdot1(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[1 + (4 * 2)] = errorFunc3.JacobianXdot2(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[2 + (4 * 2)] = errorFunc3.JacobianXdot3(meanChord, Xref, Xdotref, mass, deltaEref);
	AA[3 + (4 * 2)] = errorFunc3.JacobianXdot4(meanChord, Xref, Xdotref, mass, deltaEref);

	AA[0 + (4 * 3)] = errorFunc4.JacobianXdot1(Xref, Xdotref);
	AA[1 + (4 * 3)] = errorFunc4.JacobianXdot2(Xref, Xdotref);
	AA[2 + (4 * 3)] = errorFunc4.JacobianXdot3(Xref, Xdotref);
	AA[3 + (4 * 3)] = errorFunc4.JacobianXdot4(Xref, Xdotref);

	std::vector<double> BB;
	BB.resize(4);

	BB[0] = errorFunc1.JacobianU1(meanChord, Xref, Xdotref, mass, deltaEref);
	
	BB[1] = errorFunc2.JacobianU1(meanChord, Xref, Xdotref, mass, deltaEref);
	
	BB[2] = errorFunc3.JacobianU1(meanChord, Xref, Xdotref, mass, deltaEref);
	
	BB[3] = errorFunc4.JacobianU1(Xref, Xdotref);


	Matriz mAA(AA, 4, 4);
	Matriz mBB(BB, 4, 1);
	Matriz mEE(EE, 4, 4);

	


    return 0;
}