
#include "physics_math.h"

// PLANE //
glm::vec4 PhysicsMath::plane(glm::vec3 * p1, glm::vec3 * p2, glm::vec3 * p3)
{
	glm::vec3 d1 = glm::vec3(p2->x - p1->x, p2->y - p1->y, p2->z - p1->z);
	glm::vec3 d2 = glm::vec3(p3->x - p1->x, p3->y - p1->y, p3->z - p1->z);
	glm::vec3 n = glm::vec3(d1.y * d2.z - d1.z * d2.y, d1.z * d2.x - d1.x * d2.z, d1.x * d2.y - d1.y * d2.x);

	glm::vec4 t;
	t.w = n.x * p1->x + n.y * p1->y + n.z * p1->z;
	t.x = n.x;
	t.y = n.y;
	t.z = n.z;
	return t;
}
glm::vec4 PhysicsMath::plane(Triangle * Triangle)
{
	return plane(Triangle->t1, Triangle->t2, Triangle->t3);
}


glm::vec3 PhysicsMath::TrianglePointsNormal(glm::vec3 * p1, glm::vec3 * p2, glm::vec3  * p3)
{
	glm::vec3 d1 = glm::vec3(p2->x - p1->x, p2->y - p1->y, p2->z - p1->z);
	glm::vec3 d2 = glm::vec3(p3->x - p1->x, p3->y - p1->y, p3->z - p1->z);
	glm::vec3 n = glm::vec3(d1.y * d2.z - d1.z * d2.y, d1.z * d2.x - d1.x * d2.z, d1.x * d2.y - d1.y * d2.x);
	return n;
}
glm::vec3 PhysicsMath::TriangleNormal(Triangle * Triangle)
{
	return TrianglePointsNormal(Triangle->t1, Triangle->t2, Triangle->t3);

}

glm::vec3 PhysicsMath::getRotatedVector(glm::vec3 q, glm::vec3 p)
{
	return glm::quat(q.r, q.x, q.y, q.z) * p;
}

 float		PhysicsMath::TriangleArea3(glm::vec3 *p0, glm::vec3 *p1, glm::vec3 *p2)
{
	float d1 = (p1->y - p0->y)*(p2->z - p0->z) - (p1->z - p0->z)*(p2->y - p0->y);
	float d2 = (p1->z - p0->z)*(p2->x - p0->x) - (p1->x - p0->x)*(p2->z - p0->z);
	float d3 = (p1->x - p0->x)*(p2->y - p0->y) - (p1->y - p0->y)*(p2->x - p0->x);

	return sqrtf(d1*d1 + d2*d2 + d3*d3) / 2;
}
 bool			PhysicsMath::isInside(glm::vec3 *p, glm::vec3 *t1, glm::vec3 *t2, glm::vec3*t3)
{
	float totArea = TriangleArea3(p, t1, t2) + TriangleArea3(p, t2, t3) + TriangleArea3(p, t3, t1);
	float actArea = TriangleArea3(t1, t2, t3);
	return totArea <= actArea*1.005f;
}
 glm::vec3			PhysicsMath::getLineTriangleIntersect(glm::vec4 * p, Triangle * Triangle, glm::vec3 q1, glm::vec3 q2)
{
	glm::vec3 d = (q2 - q1);
	float t = (p->w - p->x * q1.x - p->y * q1.y - p->z * q1.z) / (p->x * d.x + p->y * d.y + p->z * d.z);
	if (t >= 0.0f && t <= 1)
	{
		glm::vec3 intPoint = q1 + (d*t);
		if (isInside(&intPoint, Triangle->t1, Triangle->t2, Triangle->t3))
			return intPoint;
	}
	return glm::vec3();
}

PhysicsMath::Intersect PhysicsMath::get2TriangleIntVec(glm::vec3 * p1, Triangle * t1, glm::vec3 * p2, Triangle * t2)
{
	Intersect i;
	glm::vec4 p = plane(&((*t2->t1) + (*p2)), &((*t2->t2) + (*p2)), &((*t2->t3) + (*p2)));
	glm::vec3 q1 = (*p1) + (*t1->t1);
	glm::vec3 q2 = (*p1) + (*t1->t2);
	glm::vec3 d = (q2 - q1);
	float t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
	if (t >= 0.0f && t <= 1)
	{
		glm::vec3 intPoint = (d*t) + q1;
		i.point = intPoint;
		if (isInside(&intPoint, t2->t1, t2->t2, t2->t3))
		{
			i.hit = 1;
			return i;
		}
	}
	q1 = *p1 + *t1->t2;
	q2 = *p1 + *t1->t3;
	d = (q2 - q1);
	t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
	if (t >= 0.0f && t <= 1)
	{
		glm::vec3 intPoint = (d*t) + q1;
		i.point = intPoint;
		if (isInside(&intPoint, t2->t1, t2->t2, t2->t3))
		{
			i.hit = 1;
			return i;
		}
	}
	q1 = *p1 + *t1->t3;
	q2 = *p1 + *t1->t1;
	d = (q2 - q1);
	t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
	if (t >= 0.0f && t <= 1)
	{
		glm::vec3 intPoint = (d*t) + q1;
		if (isInside(&intPoint, t2->t1, t2->t2, t2->t3))
			i.point = intPoint;
		if (isInside(&intPoint, t2->t1, t2->t2, t2->t3))
		{
			i.hit = 1;
			return i;
		}
	}
	return  i;
}

 float		PhysicsMath::pointPlaneDistance(glm::vec4 * plane, glm::vec3 * point)
{
	return (plane->x * point->x + plane->y * point->y + plane->z * point->z + plane->w) / plane->length();
}

 float * PhysicsMath::invert(float * m)
{
	float * i = new float[9];

	float det =
		m[0] * (m[4] * m[8] - m[5] * m[7]) -
		m[1] * (m[3] * m[8] - m[5] * m[6]) +
		m[2] * (m[3] * m[7] - m[4] * m[6]);
	i[0] = (m[4] * m[8] - m[5] * m[7]) / det;
	i[1] = -(m[1] * m[8] - m[2] * m[7]) / det;
	i[2] = (m[1] * m[5] - m[2] * m[4]) / det;
	i[3] = -(m[3] * m[8] - m[5] * m[6]) / det;
	i[4] = (m[0] * m[8] - m[2] * m[6]) / det;
	i[5] = -(m[0] * m[5] - m[2] * m[3]) / det;
	i[6] = (m[3] * m[7] - m[4] * m[6]) / det;
	i[7] = -(m[0] * m[7] - m[1] * m[6]) / det;
	i[8] = (m[0] * m[4] - m[1] * m[3]) / det;

	return i;
}
