
//// Version 1.0
//// Basic Template for a 3d OpenGL Application
#include "stdafx.h"
#include "CustomTools.h"

using namespace Tools;

#include <vector>
#include <minmax.h>

#define INF_MASS 1e32

const float g_f =  -0.001;
const float e = 0.5;
 vec4 g = vec4(0,g_f,0);
const float glob_tens = 1;
const float air_friction = 1;

 float dt = 1;

vec4 cube[8] = {vec4(), vec4(1,0,0), vec4(1,0,1), vec4(0,0,1), vec4(0,1,0), vec4(1,1,0), vec4(1,1,1), vec4(0,1,1) };
vec4 q = vec4(pi / 2, 0, 1, 0);
vec4 tetra[4] = {vec4(), vec4(1,0,0), vec4(0,0,1), vec4(0,1,0)};
  
const vec4 null_vec;//(INT_MAX, INT_MAX, INT_MAX);

vec4 g_p1;

struct point
{
	vec4 pos, vel;

	point() {};

	point(vec4 pos)
	{
		this->pos = pos;
		this->vel = vec4();
	}

	point(vec4 pos, vec4 vel)
	{
		this->pos = pos;
		this->vel = vel;
	}

	void update()
	{
		vel += g;
		vel *= air_friction;
		pos += vel;
	}
};
struct tensor
{
	point * p1, * p2;
	vec4 dir;

	float length = 0;
	float tension = 0;

	tensor() {}

	tensor(float tension, point * p1, point * p2)
	{
		this->length = ((*p2).pos - (*p1).pos).magnitude();
		this->tension = tension;
		this->p1 = p1;
		this->p2 = p2;
	}

	void calc()
	{
		vec4 dVec = (*p2).pos - (*p1).pos;
		float d = dVec.magnitude();
		dVec.normalize();
		float offset = (d - length) / 2;
		float force = offset * tension;
		dir = (dVec)* force;
		
	}
	void update()
	{
		(*p1).vel += dir;
		(*p2).vel -= dir;
	}
};
struct body
{
	bool is_static = false;
	vector<point> points;
	vector<tensor> tensors;

	body(){}

	body(vector<point> * points, vector<tensor> * tensors)
	{
		this->points = *points;
		this->tensors = *tensors;
	}

	void update()
	{
		if (!is_static)
		{
			
			for (int i = 0; i < tensors.size(); ++i)
				tensors[i].calc();
			for (int i = 0; i < tensors.size(); ++i)
				tensors[i].update();
			for (int i = 0; i < points.size(); ++i)
				points[i].update();
		}
	}

	void is_static_set(bool b)
	{
		is_static = b;
	}

	void draw()
	{
		push();
		glBegin(GL_LINES);
		for (int i = 0; i < tensors.size(); ++i)
		{
			tensors[i].p1->pos.glVertex3();
			tensors[i].p2->pos.glVertex3();
		}
		ge();
		pop();
	}
};

vec4 plane(vec4 * p1, vec4 * p2, vec4  * p3)
{
	vec4 d1 = vec4(p2->x - p1->x, p2->y - p1->y, p2->z - p1->z);
	vec4 d2 = vec4(p3->x - p1->x, p3->y - p1->y, p3->z - p1->z);
	vec4 n = vec4(d1.y * d2.z - d1.z * d2.y, d1.z * d2.x - d1.x * d2.z, d1.x * d2.y - d1.y * d2.x);

	vec4 t;
	t.w = n.x * p1->x + n.y * p1->y + n.z * p1->z;
	t.x = n.x;
	t.y = n.y;
	t.z = n.z;
	return t;
}
vec4 plane(tri * tri)
{
	return plane(tri->t1, tri->t2, tri->t3);
}

vec4 triPointsNormal(vec4 * p1, vec4 * p2, vec4  * p3)
{
	vec4 d1 = vec4(p2->x - p1->x, p2->y - p1->y, p2->z - p1->z);
	vec4 d2 = vec4(p3->x - p1->x, p3->y - p1->y, p3->z - p1->z);
	vec4 n = vec4(d1.y * d2.z - d1.z * d2.y, d1.z * d2.x - d1.x * d2.z, d1.x * d2.y - d1.y * d2.x);
	return n;
}
vec4 triNormal(tri * tri)
{
	return triPointsNormal(tri->t1, tri->t2, tri->t3);

}

vec4 getRotatedVector(vec4 q, vec4 p)
{
	return q.getRotatedPoint(p) - p;
}
inline vec4			getLineTriIntersect(vec4 * p, tri * tri, vec4 q1, vec4 q2)
{
	vec4 d = (q2 - q1);
	float t = (p->w - p->x * q1.x - p->y * q1.y - p->z * q1.z) / (p->x * d.x + p->y * d.y + p->z * d.z);
	if (t >= 0.0f && t <= 1)
	{
		vec4 intPoint = q1 + (d*t);
		if (isInside(&intPoint, tri->t1, tri->t2, tri->t3))
			return intPoint;
	}
	return vec4();
}
struct Intersect
{
	bool hit = 0;
	vec4 point;
	Intersect() {};
	Intersect(bool hit, vec4 point)
	{
		this->hit = hit;
		this->point = point;
	}
};
Intersect get2TriIntVec(vec4 * p1, tri * t1, vec4 * p2, tri * t2)
{
	Intersect i;
	vec4 p = plane(&((*t2->t1) + (*p2)), &((*t2->t2) + (*p2)), &((*t2->t3) + (*p2)));
	vec4 q1 = (*p1) + (*t1->t1);
	vec4 q2 = (*p1) + (*t1->t2);
	vec4 d = (q2 - q1);
	float t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
	if (t >= 0.0f && t <= 1)
	{
		vec4 intPoint = (d*t) + q1;
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
		vec4 intPoint = (d*t) + q1;
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
		vec4 intPoint = (d*t) + q1;
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

inline float		pointPlaneDistance(vec4 * plane, vec4 * point)
{
	return (plane->x * point->x + plane->y * point->y + plane->z * point->z + plane->w) / plane->magnitude();
}

inline float * invert(float * m)
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

inline vec4 m33v(float * m, vec4 v)
{
	return vec4(
		m[0] * v.x + m[1] * v.y + m[2] * v.z,
		m[3] * v.x + m[4] * v.y + m[5] * v.z,
		m[6] * v.x + m[7] * v.y + m[8] * v.z
	);
}

struct RigidBody
{
	vector<vec4> verts;
	vector<tri> tris;
	vec4 rotQuart, rotVec;
	vec4 pos, vel;

	float * I = new float[9];

	float mass = .1;
	float density = 1;

	float maxR;

	bool isStatic = 0;

	RigidBody() { };
	RigidBody(vector<vec4> *verts, vector<tri> *tris)
	{
		this->verts = *verts;
		this->tris = *tris;
		setMaxR();
	}
	RigidBody(vector<vec4> *verts, vector<tri> *tris, vec4 pos, vec4 vel, vec4 q)
	{
		this->verts = *verts;
		this->tris = *tris;
		this->pos = pos;
		this->vel = vel;
		this->rotVec = q;
		this->rotQuart = q.getConQuart();
		setMaxR();
	}
	RigidBody(vec4 pos, vec4 vel, vec4 q)
	{
		this->pos = pos;
		this->vel = vel;
		this->rotVec = q;
		this->rotQuart = q.getConQuart();
	}

	float * setInertiaTensor()
	{
		float * ii = new float[9];
		for (int i = 0; i < 9; ++i)
			ii[i] = 0;			
		for (int i = 0; i < verts.size(); ++i)
		{
			ii[0] += verts[i].y*verts[i].y + verts[i].z*verts[i].z;
			ii[1] += -(verts[i].x * verts[i].y);
			ii[2] += -(verts[i].x * verts[i].z);
			ii[3] += -(verts[i].x * verts[i].y);
			ii[4] += verts[i].x*verts[i].x + verts[i].z*verts[i].z;
			ii[5] += -(verts[i].y * verts[i].z);
			ii[6] += -(verts[i].x * verts[i].z);
			ii[7] += -(verts[i].y * verts[i].z);
			ii[8] += verts[i].x*verts[i].x + verts[i].y*verts[i].y;
		}
		for (int i = 0; i < 9; ++i)
			ii[i] *= mass / 12;
		return	(I = ii);
	}

	void setMaxR()
	{
		for (int i = 0; i < verts.size(); ++i)
		{
			float r = (pos - verts[i]).magnitude();
			if (r > maxR)
				maxR = r;
		}
	}

	void rotate()
	{
		vec4 r = rotVec;
		r.w *= dt;
		(rotQuart = r.getConQuart()).rotateArray(&verts);
		//rotQuart.rotateArray(&verts);
	}
	void tickGravity()
	{
		vel += g * dt;
	}
	void update()
	{
		if (!isStatic)
		{
			rotate();
			tickGravity();
			pos += vel * dt;
		}
	}

	struct VertTri
	{
		vec4 * vert;
		tri * tri;
		float d;
	};

	void checkCollision(RigidBody r)
	{
		// macro test
		if ((pos - r.pos).magnitude() > maxR + r.maxR)
		{
			//cout << "not close enough" << endl;
			return;
		}

		// micro test 

		// collision detection
		float furthest = 0;
		vec4 * vert = nullptr;
		tri * tri = nullptr;
		vector<VertTri> vertTris;
		// find vertex furthest into the body r
		// for each tri of r
		for (int j = 0; j < r.tris.size(); ++j)
		{
			vec4 p = plane(&r.tris[j]);
			// for each vert of this
			for (int i = 0; i < verts.size(); ++i)
			{
				vec4 pt = verts[i] + pos;
				// find intersect of vert to pos and the plane of the tri
				vec4 q1 = verts[i];
				vec4 q2 = pos;
				vec4 d = (pos - verts[i]);
				float t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
				// if behind the tri
				if (t >= 0.0f && t <= 1)
				{
					vec4 intPoint = q1 + (d*t);
					// if the point of intersection is inside of the tri
					if (isInside(&intPoint, (&r.tris[j])->t1, (&r.tris[j])->t2, (&r.tris[j])->t3))
					{
						t = pointPlaneDistance(&p, &pt);
					
						// update furthest point in
						if (t < furthest)
						{
							vert = &verts[i];
							tri = &r.tris[j];
							furthest = t;
							VertTri vt;
							vt.d = t;
							vt.tri = &r.tris[j];
							vt.vert = &verts[i];
							vertTris.push_back(vt);
						}
					}
				}
			}
		}

		for (int i = 0; i < vertTris.size(); ++i)
		{
			vert = vertTris[i].vert;
			tri = vertTris[i].tri;
			// collision response
			if (vert != nullptr && tri != nullptr)
			{
				//cout << "Collision!" << endl;

				// push from body r

				// normal of the tri that is bouncing off
				vec4 n = triNormal(tri).qNormalized();
				// vector from furthest point to plane orthoganol
				vec4 nt = n * furthest;
				// push
				pos -= nt;
				// collision point on tri/plane
				vec4 p = (*vert + pos);



				// apply impulse
				float mass = 1;
				float mu = 0.5;

				vec4 r1 = p - r.pos;
				vec4 r2 = *vert;

				vec4 w1 = getRotatedVector(r.rotQuart, p - r.pos);
				vec4 w2 = getRotatedVector(rotQuart, p - pos);

				vec4 i1 = w1.normalized();
				vec4 i2 = w2.normalized();

				vec4 vr = vel - r.vel;

				float top = vr.dot(n) * -(1 + e);
				vec4 left = (r1.Cross(n).Cross(r1) / i1);
				left.normalize();
				vec4 right = (r2.Cross(n).Cross(r2) / i2);
				right.normalize();
				float bottom = (left + right).dot(n);
				float jr = max(top / bottom,0);

				vec4 impulse1 = n * -jr;
				vec4 impulse2 = n * jr;

				//r.vel += impulse1.inverse();
				vel += impulse2;
				

				vec4 W1 = w1 - ((r1.Cross(n) / i1) * jr);
				vec4 W2 = w2 + ((r2.Cross(n) / i2) * jr);

				vec4 crossPVelNorm = W2.Cross(r2);
				rotVec = crossPVelNorm.normalized();
				rotVec.w = W2.magnitude();

				//vec4 D = pos - *vert;
				//vec4 u = getRotVec(rotQuart, *vert);
				//vec4 S = vel + u;

				//vec4 R = normal * (S.magnitude() * S.Dot(normal * -1));
				//vec4 F = R * mu;

				//F = S.Cross(D).Cross(S).normalized() * F.magnitude();

				//float j = max(-(1 + e) * S.Dot(normal), 0);
				//vec4 SP = S + (normal * (S.magnitude() * S.Dot(normal) * e));				
				//SP -= F;

				//vec4 T = D.normalized() * (SP.magnitude() * D.Dot(SP));
				//vel -= T;
				////vel -= F;

				////SP -= F;
				//// W needs to be projected onto the perp to D but int the direction of V
				//vec4 W = SP.Cross(D).Cross(SP).normalized() * -SP.magnitude();// v.Cross(D).Cross(v).normalized() * v.magnitude();
				//
				//vec4 crossPVelNorm = W.Cross(D).normalized();
				//rotVec = crossPVelNorm.normalized();
				//rotVec.w = W.magnitude();

				//vec4 vertD = (pos - *vert).magnitude();
				//float bodyImpulseRatio = normal.Dot(vertD);

				//vec4 oldV = vel + getRotVec(rotQuart, *vert);
				//float reflectRatio = normal.Dot(oldV * -1);
				//float j = max(-(1 + e) * reflectRatio, 0);
				//vec4 newV = oldV * -e;// + normal * (j * oldV.magnitude());

				////j = max(-(1 + e) * bodyImpulseRatio, 0);
				//vel += normal * (reflectRatio * oldV.magnitude());






				//// apply impulse
				//float mass = 1;
				//float d = vel.normalized().dot(normal);
				//float j = max(-(1 + e) * d, 0);
				////vel += normal * (j * vel.magnitude());

				//vec4 vertD = (pos - *vert).magnitude();
				//float t = normal.dot(vertD);
				//vec4 oldV = vel + getRotVec(rotQuart, *vert);
				//d = oldV.dot(normal);
				////j = max(-(1 + e) * d, 0);
				//vec4 newV = normal * (d * oldV.magnitude());

				//vec4 crossPVelNorm = newV.Cross(vertD).normalized();
				//rotVec = crossPVelNorm.normalized();
				//rotVec.w = newV.magnitude();

				//vel += normal * (t * newV.magnitude());


			}
		}
	}

	void checkCollision2(RigidBody r)
	{
		bool finish = 0;
		//for (int i = 0; i < tris.size(); ++i)
		//{
		//	for (int j = 0; j < r.tris.size(); ++j)
		//	{
		//		Intersect t = get2TriIntVec(&pos, &tris[i], &r.pos, &r.tris[j]);
		//		if (t.hit)
		//		{
		//			
		//			pos += ;
		
		//			////cout << t.point.toString() << endl;
		//			//t.point.glTranslate();
		//			//if(0)
		//			//if (r.isStatic)
		//			//{
		//			//	rotVec *= -1;
		//			//	//cout << "hit" << endl;;
		//			//}
		//			//else
		//			//{
		//			//	vec4 v1 = rotVec.Cross(pos - t.point);
		//			//	vec4 v2 = r.rotVec.Cross(r.pos - t.point);
		//			//	vec4 newV = v1.Cross(v2);
		//			//	newV.normalize();
		//			//	rotVec += newV;
		//			//	rotVec += getRotVec(rotQuart, t.point - pos).Cross(getRotVec(r.rotQuart, t.point - r.pos));
		//			//	rotVec.normalize();
		//			//}
		//			////rotate();
		//			//float d = (pos - t.point).magnitude();
		//			////vec4 v = rotVec.Cross(pos - t.point);
		//			////rotVec.w *= d;
		//			////rotate();
		//			////rotQuart = rotVec.getConQuart();
		//			////vel += (getRotVec(rotQuart, t.point - pos)) / d;
		//			//vel += vel * -2 /d;
		//			//pos += vel;
		//			//finish = 1;
		//		}
		//		if (finish) break;
		//	}
		//	if (finish) break;
		//}

		// collision detection
		float furthest = 0;
		vec4 * vert = nullptr;
		tri * tri = nullptr;
		// find vertex furthest into the body r
		for (int j = 0; j < r.tris.size(); ++j)
		{
			vec4 p = plane(&r.tris[j]);
			for (int i = 0; i < verts.size(); ++i)
			{
				vec4 pt = verts[i] + pos;
				if (getLineTriIntersect(&p, &r.tris[j],verts[i],pos) != null_vec)
				{
					float t = pointPlaneDistance(&p, &pt);
					if (t < furthest)
					{
						vert = &verts[i];
						tri = &r.tris[j];
						furthest = t;
					}
				}
			}				
		}
		// collision response
		if (vert != nullptr && tri != nullptr)
		{
			cout << "Collision!!" << endl;

			// push from body r

			// normal of the tri that is bouncing off
			vec4 normal = triNormal(tri).qNormalized();
			// vector from furthest point to plane orthoganol
			vec4 nt = normal * furthest;
			// push
			pos -= nt;
			// collision point on tri/plane
			vec4 colPoint = (*vert + pos) + nt;

			// apply impulse
			float mass = 1;
			float d = vel.normalized().dot(normal);
			float j = max(-(1 + e) * d, 0);

			//cout << j << endl;
			vec4 cVel = vel;// +getRotVec(rotQuart, *vert);
			vel += normal * (j * vel.magnitude());


			/*
			vec4 crossPVelNorm = cVel.Cross(pos - *vert);
			rotVec = crossPVelNorm;
			float r = (pos - *vert).magnitude();
			rotVec.w -= cVel.magnitude() / r / 5;
			*/

			//vec4 crossPVelNorm = cVel.Cross(pos - *vert);
			//rotVec = crossPVelNorm;
			//float r = (pos - *vert).magnitude();
			//rotVec.w -= cVel.magnitude() / r / 5;

		}

	}

	void drawVerts()
	{
		push();
		pos.glTranslate();
		glBegin(GL_LINE_LOOP);
		for (int i = 0; i < tris.size(); ++i)
			tris[i].draw();
		ge();
		pop();
	}
	void drawTris()
	{
		push();
		pos.glTranslate();
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < tris.size(); ++i)
		{
			float r = (int(255 * i / (float)tris.size()) % 255) / 255.0f;
			float g = (int(255 * (i + tris.size() / 3) / (float)tris.size()) % 255) / 255.0f;
			float b = (int(255 * (i + tris.size() * 2 / 3) / (float)tris.size()) % 255) / 255.0f;
			vec4(r,g,b).setAlpha(1).glColour4();
			tris[i].draw();
		}
		ge();
		pop();
	}
	void drawNormals()
	{
		COLOUR_GREEN.glColour4();
		push();
		pos.glTranslate();
		glBegin(GL_LINES);
		for (int i = 0; i < tris.size(); ++i)
		{
			vec4 n = triNormal(&tris[i]);
			vec4 p = (*tris[i].t1 + (*tris[i].t2 + *tris[i].t3));
			p /= 3;
			//p = (p + *tris[i].t1) / 2;
			//p = *tris[i].t2;
			p.glVertex3();
			(p + (n*2)).glVertex3();
		}
		ge();
		pop();
	}
};

void checkCollision(RigidBody* b1, RigidBody* b2)
{
	// macro test
	if ((b1->pos - b2->pos).magnitude() > b1->maxR + b2->maxR)
	{
		//cout << "not close enough" << endl;
		return;
	}

	// micro test 

	// collision detection
	float furthest = 0;
	vec4 * vert = nullptr;
	tri * tri = nullptr;
	//vector<VertTri> vertTris;
	// find vertex furthest into the body r
	// for each tri of r
	for (int j = 0; j < b2->tris.size(); ++j)
	{
		vec4 p = plane(&b2->tris[j]);
		// for each vert of this
		for (int i = 0; i < b1->verts.size(); ++i)
		{
			vec4 pt = b1->verts[i] + b1->pos;
			// find intersect of vert to pos and the plane of the tri
			vec4 q1 = b1->verts[i];
			vec4 q2 = b1->pos;
			vec4 d = (b1->pos - b1->verts[i]);
			float t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
			// if behind the tri
			if (t >= 0.0f && t <= 1)
			{
				vec4 intPoint = q1 + (d*t);
				// if the point of intersection is inside of the tri
				if (isInside(&intPoint, (&b2->tris[j])->t1, (&b2->tris[j])->t2, (&b2->tris[j])->t3))
				{
					t = pointPlaneDistance(&p, &pt);

					// update furthest point in
					if (t < furthest)
					{
						vert = &b1->verts[i];
						tri = &b2->tris[j];
						furthest = t;
						//VertTri vt;
						//vt.d = t;
						//vt.tri = &r.tris[j];
						//vt.vert = &verts[i];
						//vertTris.push_back(vt);
					}
				}
			}
		}
	}

	//for (int i = 0; i < vertTris.size(); ++i)
	{
		//vert = vertTris[i].vert;
		//tri = vertTris[i].tri;

		// collision response
		if (vert != nullptr && tri != nullptr)
		{
			//cout << "Collision!" << endl;

			// push from body r

			// normal of the tri that is bouncing off
			vec4 n = triNormal(tri).normalized();
			// vector from furthest point to plane orthoganol
			vec4 nt = n * furthest;
			// push	
			if (!b1->isStatic && !b2->isStatic)
			{
				//cout << 1 << endl;
				b1->pos -= nt;
				//b2->pos += nt / 2;
			}
			else if (b1->isStatic)
			{
				//cout << 2 << endl;

				//b2->pos += nt;
			}
			else if (b2->isStatic)
			{
				//cout << 3 << endl;

				b1->pos -= nt;
			}

			// collision point on tri/plane
			vec4 p = (*vert + b1->pos);

			//cout << b1->vel.magnitude() << endl;

			// calc impulse
			float mu = 1;

			// collision point vectors
			vec4 r1 = p - b1->pos;
			vec4 r2 = p - b2->pos;

			// inertia tensors
			float * i1_1 = invert(b1->I); //r1 * b1->rotVec.w;//w1.normalized();
			//float * i2_1 = invert(b2->I); // r2 * b2->rotVec.w;//w2.normalized();
			
			float i1 = vec4(
				n.x * b1->I[0] + n.y * b1->I[3] + n.z * b1->I[6],
				n.x * b1->I[1] + n.y * b1->I[4] + n.z * b1->I[7],
				n.x * b1->I[2] + n.y * b1->I[5] + n.z * b1->I[8]
			).dot(n);
			float i2 = vec4(
				n.x * b2->I[0] + n.y * b2->I[3] + n.z * b2->I[6],
				n.x * b2->I[1] + n.y * b2->I[4] + n.z * b2->I[7],
				n.x * b2->I[2] + n.y * b2->I[5] + n.z * b2->I[8]
			).dot(n);

			// relative vel of CoMs
			vec4 vr = b2->vel - b1->vel;

			// impulse ratio

			float jr = 
				(-(1+e) * vr.dot(n))
				/ 
				(1/b1->mass + 1/b2->mass + 
					(r1.Cross(n).dot(r1.Cross(n))/i1) + 
					(r2.Cross(n).dot(r2.Cross(n)) / i2));
			jr *= mu;

			// impulse vectors
			vec4 impulse1 = n * (-jr / b1->mass);
			//vec4 impulse2 = n * (jr / b2->mass);

			// apply impulse
			b1->vel += impulse1;
			//b2->vel += impulse2;

			// angular vector of collision point
			vec4 w1 = getRotatedVector(b1->rotQuart, r1);
			//vec4 w2 = getRotatedVector(b2->rotQuart, r2);

			// angular impulse vectors
			vec4 W1 = w1 - (r1.Cross(n)* (jr/i1));

			//vec4 W2 = w2 + (m33v(i2_1, r2.Cross(n)) * jr);

			// calc new quarternion for rotation and 
			// apply impulse magnitude to the angular speed
			//b1->rotVec = W1;
			//b1->rotVec.w = W1.magnitude() / r1.magnitude();
			//b2->rotVec = W2.Cross(r2).normalized();
			//b2->rotVec.w = W2.magnitude() / r1.magnitude();

		}
	}
}


void checkCollision3(RigidBody* b1, RigidBody* b2)
{
	// macro test
	if ((b1->pos - b2->pos).magnitude() > b1->maxR + b2->maxR)
	{
		//cout << "not close enough" << endl;
		return;
	}

	// micro test 

	// collision detection
	float furthest = 0;
	vec4 * vert = nullptr;
	tri * tri = nullptr;
	//vector<VertTri> vertTris;
	// find vertex furthest into the body r
	// for each tri of r
	for (int j = 0; j < b2->tris.size(); ++j)
	{
		vec4 p = plane(&b2->tris[j]);
		// for each vert of this
		for (int i = 0; i < b1->verts.size(); ++i)
		{
			vec4 pt = b1->verts[i] + b1->pos;
			// find intersect of vert to pos and the plane of the tri
			vec4 q1 = b1->verts[i];
			vec4 q2 = b1->pos;
			vec4 d = (b1->pos - b1->verts[i]);
			float t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
			// if behind the tri
			if (t >= 0.0f && t <= 1)
			{
				vec4 intPoint = q1 + (d*t);
				// if the point of intersection is inside of the tri
				if (isInside(&intPoint, (&b2->tris[j])->t1, (&b2->tris[j])->t2, (&b2->tris[j])->t3))
				{
					t = pointPlaneDistance(&p, &pt);

					// update furthest point in
					if (t < furthest)
					{
						vert = &b1->verts[i];
						tri = &b2->tris[j];
						furthest = t;
						//VertTri vt;
						//vt.d = t;
						//vt.tri = &r.tris[j];
						//vt.vert = &verts[i];
						//vertTris.push_back(vt);
					}
				}
			}
		}
	}

	//for (int i = 0; i < vertTris.size(); ++i)
	{
		//vert = vertTris[i].vert;
		//tri = vertTris[i].tri;

		// collision response
		if (vert != nullptr && tri != nullptr)
		{
			//cout << "Collision!" << endl;

			// push from body r

			// normal of the tri that is bouncing off
			vec4 n = triNormal(tri).normalized();
			// vector from furthest point to plane orthoganol
			vec4 nt = n * furthest;
			// push	
			if (!b1->isStatic && !b2->isStatic)
			{
				//cout << 1 << endl;
				b1->pos -= nt / 2;
				b2->pos += nt / 2;
			}
			else if (b1->isStatic)
			{
				//cout << 2 << endl;

				b2->pos += nt;
			}
			else if (b2->isStatic)
			{
				//cout << 3 << endl;

				b1->pos -= nt;
			}

			// collision point on tri/plane
			vec4 p = (*vert + b1->pos);

			//cout << b1->vel.magnitude() << endl;

			// calc impulse
			float mu = 0.5;

			// collision point vectors
			vec4 r1 = p - b1->pos;
			vec4 r2 = p - b2->pos;

			// inertia tensors
			float * i1_1 = invert(b1->I); //r1 * b1->rotVec.w;//w1.normalized();
			float * i2_1 = invert(b2->I); // r2 * b2->rotVec.w;//w2.normalized();
									   // difference in vel
			vec4 vr = b2->vel - b1->vel;

			vec4 tv = vr - (n * vr.dot(n));
			tv.normalize();
			//vec4 R = n * (vr.dot(n) * vr.magnitude());
			//tv = tv * (R.magnitude() * mu);

			// impulse ratio
			float top = -(1 + e) * vr.dot(n) * dt;

			vec4 left = m33v(i1_1,r1.Cross(n)).Cross(r1);
			left.normalize();
			if (b1->isStatic)
				left = vec4();

			vec4 right = m33v(i2_1,r2.Cross(n)).Cross(r2);
			right.normalize();
			if (b2->isStatic)
				right = vec4();

			float bottom = (left + right).dot(n);
			bottom += 1 / b1->mass;
			bottom += 1 / b2->mass;

			float jr = top / bottom;

			// impulse vectors
			vec4 impulse1 = n * (-jr / b1->mass);
			vec4 impulse2 = n * (jr / b2->mass);

			// apply impulse
			b1->vel += impulse1;
			//b2->vel += impulse2;

			// angular vector of collision point
			vec4 w1 = getRotatedVector(b1->rotQuart, r1);
			//vec4 w2 = getRotatedVector(b2->rotQuart, r2);
			
			// angular impulse vectors
			vec4 W1 = w1 - (m33v(i1_1, r1.Cross(n*jr)));
			//vec4 W2 = w2 + (m33v(i2_1, r2.Cross(n)) * jr);

			// calc new quarternion for rotation and 
			// apply impulse magnitude to the angular speed
			b1->rotVec = W1.Cross(r1).normalized();
			b1->rotVec.w = W1.magnitude() / (r1.magnitude()*pi * 2);
			//b2->rotVec = W2.Cross(r2).normalized();
			//b2->rotVec.w = W2.magnitude() / r1.magnitude();

		}
	}
}

void checkCollision2(RigidBody* b1, RigidBody* b2)
{
	// macro test
	if ((b1->pos - b2->pos).magnitude() > b1->maxR + b2->maxR)
	{
		//cout << "not close enough" << endl;
		return;
	}

	// micro test 

	// collision detection
	float furthest = 0;
	vec4 * vert = nullptr;
	tri * tri = nullptr;
	//vector<VertTri> vertTris;
	// find vertex furthest into the body r
	// for each tri of r
	for (int j = 0; j < b2->tris.size(); ++j)
	{
		vec4 p = plane(&b2->tris[j]);
		// for each vert of this
		for (int i = 0; i < b1->verts.size(); ++i)
		{
			vec4 pt = b1->verts[i] + b1->pos;
			// find intersect of vert to pos and the plane of the tri
			vec4 q1 = b1->verts[i];
			vec4 q2 = b1->pos;
			vec4 d = (b1->pos - b1->verts[i]);
			float t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
			// if behind the tri
			if (t >= 0.0f && t <= 1)
			{
				vec4 intPoint = q1 + (d*t);
				// if the point of intersection is inside of the tri
				if (isInside(&intPoint, (&b2->tris[j])->t1, (&b2->tris[j])->t2, (&b2->tris[j])->t3))
				{
					t = pointPlaneDistance(&p, &pt);

					// update furthest point in
					if (t < furthest)
					{
						vert = &b1->verts[i];
						tri = &b2->tris[j];
						furthest = t;
						//VertTri vt;
						//vt.d = t;
						//vt.tri = &r.tris[j];
						//vt.vert = &verts[i];
						//vertTris.push_back(vt);
					}
				}
			}
		}
	}

	//for (int i = 0; i < vertTris.size(); ++i)
	{
		//vert = vertTris[i].vert;
		//tri = vertTris[i].tri;

		// collision response
		if (vert != nullptr && tri != nullptr)
		{
			//cout << "Collision!" << endl;

			// push from body r

			// normal of the tri that is bouncing off
			vec4 n = triNormal(tri).normalized();
			// vector from furthest point to plane orthoganol
			vec4 nt = n * furthest;
			// push	
			if (!b1->isStatic && !b2->isStatic)
			{
				//cout << 1 << endl;
				b1->pos -= nt / 2;
				b2->pos += nt / 2;
			}
			else if (b1->isStatic)
			{
				//cout << 2 << endl;

				b2->pos += nt;
			}
			else if (b2->isStatic)
			{
				//cout << 3 << endl;

				b1->pos -= nt;
			}

			// collision point on tri/plane
			vec4 p = (*vert + b1->pos);

			//cout << b1->vel.magnitude() << endl;

			// calc impulse
			float mu = 0.5;

			// collision point vectors
			vec4 r1 = p - b1->pos;
			vec4 r2 = p - b2->pos;

			// angular vector of collision point
			vec4 w1 = getRotatedVector(b1->rotQuart, r1);
			vec4 w2 = getRotatedVector(b2->rotQuart, r2);

			// rotation vector
			vec4 i1 = w1.normalized(); //r1 * b1->rotVec.w;//w1.normalized();
			vec4 i2 = w2.normalized(); // r2 * b2->rotVec.w;//w2.normalized();

			// difference in vel
			vec4 vr = b2->vel - b1->vel;

			vec4 tv = vr - (n * vr.dot(n));
			tv.normalize();
			//vec4 R = n * (vr.dot(n) * vr.magnitude());
			//tv = tv * (R.magnitude() * mu);

			// impulse ratio
			float top = -(1+e) * vr.dot(n) * dt;

			vec4 left = (r1.Cross(n).Cross(r1) / i1);
			left.normalize();
			if (b1->isStatic)
				left = vec4();

			vec4 right = (r2.Cross(n).Cross(r2) / i2);
			right.normalize();
			if (b2->isStatic)
				right = vec4();

			float bottom = (left + right).dot(n);
			//if (!b1->isStatic)
				bottom += 1 / b1->mass;
			//if (!b2->isStatic)
				bottom += 1 / b2->mass;
			float jr = top / bottom;

			// impulse vectors
			vec4 impulse1 = n * (-jr / b1->mass);
			//impulse1 += tv;
			vec4 impulse2 = n * (jr / b2->mass);
			//impulse1 += tv;

			// apply impulse
			b1->vel += impulse1;
			b2->vel += impulse2;

			// angular impulse vectors
			vec4 W1 = w1 - ((r1.Cross(n) / i1) * jr);
			vec4 W2 = w2 + ((r2.Cross(n) / i2) * jr);

			// calc new quarternion for rotation and 
			// apply impulse magnitude to the angular speed
			//b1->rotVec = W1.Cross(r1).normalized();
			//b1->rotVec.w = W1.magnitude()/ r1.magnitude();
			//b2->rotVec = W2.Cross(r2).normalized();
			//b2->rotVec.w = W2.magnitude() / r2.magnitude();

		}
	}
}

RigidBody rb1, rb2, ground;

body b1, b2;

vector<vec4> tracer;
Timer tracerTick(1,1);

inline void draw3D()
{
	start3dDraw();

	COLOUR_GREEN.glColour4();
	//b1.draw();
	COLOUR_RED.glColour4();
	//b2.draw();

	rb1.drawTris();
	rb1.drawNormals();
	rb2.drawTris();
	rb2.drawNormals();

	COLOUR_GREEN.glColour4();
	ground.drawTris();
	ground.drawNormals();

	COLOUR_WHITE.glColour4();

	//push();
	////glBegin(GL_POINT);
	////rb2.checkCollision(rb1);
	//glutWireCube(0.1);
	//pop();

	
	COLOUR_GREEN.glColour4();
	push();
	glBegin(GL_LINES);
	g_p1.glVertex3();
	vec4(g_p1.x, 0, g_p1.z).glTranslate();
	ge();
	pop();

	push();
	glBegin(GL_LINES);
	for (int i = 1; i < tracer.size(); ++i)
	{
		tracer[i - 1].glVertex3();
		tracer[i].glVertex3();
	}
	ge();
	pop();

	/*push();
	glBegin(GL_LINE_LOOP);
	for (int i = 0; i < 8; ++i)
		cube[i].glVertex3();
	ge();
	pop();*/

	//push();
	//poly();
	//vec4().glVertex3();
	//vec4(5, 0).glVertex3();
	//vec4(5, 0, 5).glVertex3();
	//vec4(0, 0, 5).glVertex3();
	//ge();
	//pop();

}
inline void drawDebugInfo()
{
	COLOUR_WHITE.glColour3();
	drawFPS();
}
inline void draw2D()
{
	start2dDraw();

	//gui->draw();
	drawDebugInfo();

	end2dDraw();
}

void renderScene()
{
	std::this_thread::sleep_for(std::chrono::milliseconds(16));

	glClearColor(0.1, 0.1, 0.1, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	draw3D();

	draw2D();

	glFlush();
	glutSwapBuffers();
	glutPostRedisplay();

	incFFC();
}
void processThread()
{
	while (1)
	{
		startProcessLoop();

		//for (int i = 0; i < b2.points.size(); ++i)
		//{
		//	vec4 p = b2.points[i].pos;
		//	vec4 v = b2.points[i].vel;
		//	if (p.y < 0)
		//	{
		//		v = v * -0.01;
		//		p.y = 0;
		//	}
		//	b2.points[i].pos = p;
		//	b2.points[i].vel = v;
		//}


		//b1.update();
		//b2.update();

		//worldPosition = b2.points[0].pos + vec4(0, 1, -5);
		//worldPosition.y = 1;

		//vec4(0, 1, 0, pi / 100).conQuart().rotateArray(cube,8);

		checkCollision(&rb1, &ground);		
		//checkCollision(&rb2, &ground);
		//checkCollision(&rb1, &rb2);
		//checkCollision(&rb2, &rb1);
		rb1.update();
		//rb2.update();
		//rb1.checkCollision(ground);
		//rb2.update();
		//ground.update();

		if (tracerTick.tick())
			tracer.push_back(rb1.pos);

		// cout << get2TriIntVec(&rb1.pos, &rb1.tris[0], &ground.pos, &ground.tris[0]).toString() << endl;

		endProcessLoop();
	}
}

void reshape3D(int w, int h)
{
	window = vec4(w, h);
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(fov, (float)w / (float)h, 0.1, 1000);
}
void customGLInit()
{
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(200, 100);
	glutInitWindowSize(window.x, window.y);
	glutCreateWindow(windowName.c_str());
	glutDisplayFunc(renderScene);
	glutReshapeFunc(reshape3D);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glShadeModel(GL_SMOOTH);
}

void setup_bodies()
{
	//b1.points.push_back(point(vec4()));
	//b1.points.push_back(point(vec4(0, 0, 1)));
	//b1.points.push_back(point(vec4(1, 0, 1)));
	//b1.points.push_back(point(vec4(1, 0, 0)));

	//b1.tensors.push_back(tensor(glob_tens, &b1.points[0], &b1.points[1]));
	//b1.tensors.push_back(tensor(glob_tens, &b1.points[1], &b1.points[2]));
	//b1.tensors.push_back(tensor(glob_tens, &b1.points[2], &b1.points[3]));
	//b1.tensors.push_back(tensor(glob_tens, &b1.points[3], &b1.points[0]));

	//b1.is_static_set(true);

	//b2.points.push_back(point(vec4(0, 0, 0)));
	//b2.points.push_back(point(vec4(0, 0, 1)));
	//b2.points.push_back(point(vec4(1, 0, 1)));
	//b2.points.push_back(point(vec4(1, 0, 0)));
	//b2.points.push_back(point(vec4(0, 1, 0)));
	//b2.points.push_back(point(vec4(0, 1, 1)));
	//b2.points.push_back(point(vec4(1, 1, 1)));
	//b2.points.push_back(point(vec4(1, 1, 0)));

	//b2.tensors.push_back(tensor(glob_tens, &b2.points[0], &b2.points[1]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[1], &b2.points[2]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[2], &b2.points[3]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[3], &b2.points[0]));

	//b2.tensors.push_back(tensor(glob_tens, &b2.points[4], &b2.points[5]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[5], &b2.points[6]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[6], &b2.points[7]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[7], &b2.points[4]));

	//b2.tensors.push_back(tensor(glob_tens, &b2.points[0], &b2.points[4]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[1], &b2.points[5]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[2], &b2.points[6]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[3], &b2.points[7]));

	//b2.tensors.push_back(tensor(glob_tens, &b2.points[0], &b2.points[6]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[1], &b2.points[7]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[2], &b2.points[4]));
	//b2.tensors.push_back(tensor(glob_tens, &b2.points[3], &b2.points[5]));


	//for (int i = 0; i < 4; ++i)
	//	cout << &ts[i].p1 << " : " << &b2.points[i] << endl;

	/*vector<vec4> verts;
	verts.push_back(vec4());
	verts.push_back(vec4(1,0,0));
	verts.push_back(vec4(0,0,1));
	verts.push_back(vec4(0,1,0));
*/


	rb1 = RigidBody(vec4(0,5,0), vec4(0,0,0), vec4(randf(), randf(), randf(), 0.01));
	
	rb1.verts.push_back(vec4(-1,-1,-1));
	rb1.verts.push_back(vec4(-1,-1,1));
	rb1.verts.push_back(vec4(1,-1,1));
	rb1.verts.push_back(vec4(1,-1,-1));

	rb1.verts.push_back(vec4(-1, 1, -1));
	rb1.verts.push_back(vec4(-1, 1, 1));
	rb1.verts.push_back(vec4(1, 1, 1));
	rb1.verts.push_back(vec4(1, 1, -1));

	rb1.tris.push_back(tri(&rb1.verts[0], &rb1.verts[3], &rb1.verts[2]));
	rb1.tris.push_back(tri(&rb1.verts[2], &rb1.verts[1], &rb1.verts[0]));
	rb1.tris.push_back(tri(&rb1.verts[4], &rb1.verts[5], &rb1.verts[6]));
	rb1.tris.push_back(tri(&rb1.verts[6], &rb1.verts[7], &rb1.verts[4]));

	rb1.tris.push_back(tri(&rb1.verts[0], &rb1.verts[1], &rb1.verts[5]));
	rb1.tris.push_back(tri(&rb1.verts[5], &rb1.verts[4], &rb1.verts[0]));
	rb1.tris.push_back(tri(&rb1.verts[2], &rb1.verts[3], &rb1.verts[7]));
	rb1.tris.push_back(tri(&rb1.verts[7], &rb1.verts[6], &rb1.verts[2]));
	
	rb1.tris.push_back(tri(&rb1.verts[0], &rb1.verts[4], &rb1.verts[7]));
	rb1.tris.push_back(tri(&rb1.verts[7], &rb1.verts[3], &rb1.verts[0]));
	rb1.tris.push_back(tri(&rb1.verts[1], &rb1.verts[2], &rb1.verts[6]));
	rb1.tris.push_back(tri(&rb1.verts[6], &rb1.verts[5], &rb1.verts[1]));

	rb1.setInertiaTensor();
	rb1.setMaxR();

	rb2 = RigidBody(vec4(0, 10, 0), vec4(), vec4(randf(), randf(), randf(), 0.01));
	rb2.isStatic = 1;
	rb2.verts.push_back(vec4(-1, -1, -1));
	rb2.verts.push_back(vec4(-1, -1, 1));
	rb2.verts.push_back(vec4(1, -1, 1));
	rb2.verts.push_back(vec4(1, -1, -1));

	rb2.verts.push_back(vec4(-1, 1, -1));
	rb2.verts.push_back(vec4(-1, 1, 1));
	rb2.verts.push_back(vec4(1, 1, 1));
	rb2.verts.push_back(vec4(1, 1, -1));

	rb2.tris.push_back(tri(&rb2.verts[0], &rb2.verts[3], &rb2.verts[2]));
	rb2.tris.push_back(tri(&rb2.verts[2], &rb2.verts[1], &rb2.verts[0]));
	rb2.tris.push_back(tri(&rb2.verts[4], &rb2.verts[5], &rb2.verts[6]));
	rb2.tris.push_back(tri(&rb2.verts[6], &rb2.verts[7], &rb2.verts[4]));

	rb2.tris.push_back(tri(&rb2.verts[0], &rb2.verts[1], &rb2.verts[5]));
	rb2.tris.push_back(tri(&rb2.verts[5], &rb2.verts[4], &rb2.verts[0]));
	rb2.tris.push_back(tri(&rb2.verts[2], &rb2.verts[3], &rb2.verts[7]));
	rb2.tris.push_back(tri(&rb2.verts[7], &rb2.verts[6], &rb2.verts[2]));

	rb2.tris.push_back(tri(&rb2.verts[0], &rb2.verts[4], &rb2.verts[7]));
	rb2.tris.push_back(tri(&rb2.verts[7], &rb2.verts[3], &rb2.verts[0]));
	rb2.tris.push_back(tri(&rb2.verts[1], &rb2.verts[2], &rb2.verts[6]));
	rb2.tris.push_back(tri(&rb2.verts[6], &rb2.verts[5], &rb2.verts[1]));

	rb2.setInertiaTensor();
	rb2.setMaxR();

	ground = RigidBody(vec4(), vec4(), vec4(1, 1, 1, 0));
	//ground = RigidBody(vec4(),vec4(), vec4(1));
	ground.isStatic = 1;
	ground.verts.push_back(vec4(-10, 0, -10));
	ground.verts.push_back(vec4(-10, 0, 10));
	ground.verts.push_back(vec4(10, 0, 10));
	ground.verts.push_back(vec4(10, 0, -10));

	ground.tris.push_back(tri(&ground.verts[0], &ground.verts[1], &ground.verts[2]));
	ground.tris.push_back(tri(&ground.verts[2], &ground.verts[3], &ground.verts[0]));

	ground.setInertiaTensor();
	ground.setMaxR();

	ground.mass = INF_MASS;
}

int main(int argc, char** argv)
{
	srand(time(0));
	glutInit(&argc, argv);
	customGLInit();
	win = FindWindowA(NULL, windowName.c_str());

	worldPosition = vec4(1, 0.5, 2) * 5;
	lookDirection = worldPosition * -1;

	setup_bodies();

	std::thread process(processThread);
	process.detach();

	detachFPSThread();

	//detachEventHandler();

	glutMainLoop();
	getc(stdin);
	return 0;
}

