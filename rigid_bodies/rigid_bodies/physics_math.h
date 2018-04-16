#pragma once

#include "glm.h"


namespace PhysicsMath {
	// Triangle //
	struct Triangle
	{
		glm::vec3 *t1, *t2, *t3;
		Triangle() {}
		Triangle(glm::vec3* t1, glm::vec3* t2, glm::vec3* t3)
		{
			this->t1 = t1;
			this->t2 = t2;
			this->t3 = t3;
		}
	}; 
	
	//  INTERSECT //
	struct Intersect
	{
		bool hit = 0;
		glm::vec3 point;
		Intersect() {};
		Intersect(bool hit, glm::vec3 point)
		{
			this->hit = hit;
			this->point = point;
		}
	};

	struct VertTri
	{
		glm::vec3 * vert;
		Triangle * tri;
		float d;
	};

	// PLANE //
	glm::vec4 plane(glm::vec3 * p1, glm::vec3 * p2, glm::vec3 * p3);
	glm::vec4 plane(Triangle * Triangle);


	glm::vec3 TrianglePointsNormal(glm::vec3 * p1, glm::vec3 * p2, glm::vec3  * p3);
	glm::vec3 TriangleNormal(Triangle * Triangle);

	glm::vec3 getRotatedVector(glm::vec3 q, glm::vec3 p);

	 float		TriangleArea3(glm::vec3 *p0, glm::vec3 *p1, glm::vec3 *p2);
	 bool			isInside(glm::vec3 *p, glm::vec3 *t1, glm::vec3 *t2, glm::vec3*t3);
	 glm::vec3			getLineTriangleIntersect(glm::vec4 * p, Triangle * Triangle, glm::vec3 q1, glm::vec3 q2);
	
	Intersect get2TriangleIntVec(glm::vec3 * p1, Triangle * t1, glm::vec3 * p2, Triangle * t2);

	 float		pointPlaneDistance(glm::vec4 * plane, glm::vec3 * point);

	 float * invert(float * m);

}