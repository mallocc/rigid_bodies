#include "rigid_body.h"

void checkCollision(RigidBody* b1, RigidBody* b2)
{
	// macro test
	if ((b1->pos - b2->pos).length() > b1->maxRadius + b2->maxRadius)
	{
		//printf("not close enough\n");
		return;
	}

	// micro test 

	// collision detection
	float furthest = 0;
	glm::vec3 * vert = nullptr;
	PhysicsMath::Triangle * tri = nullptr;
	//vector<VertPhysicsMath::Triangle> vertPhysicsMath::Triangles;
	// find vertex furthest into the body r
	// for each PhysicsMath::Triangle of r
	for (int j = 0; j < b2->tris.size(); ++j)
	{
		glm::vec4 p = plane(&b2->tris[j]);
		// for each vert of this
		for (int i = 0; i < b1->verts.size(); ++i)
		{
			glm::vec3 pt = b1->verts[i] + b1->pos;
			// find intersect of vert to pos and the plane of the PhysicsMath::Triangle
			glm::vec3 q1 = b1->verts[i];
			glm::vec3 q2 = b1->pos;
			glm::vec3 d = (b1->pos - b1->verts[i]);
			float t = (p.w - p.x * q1.x - p.y * q1.y - p.z * q1.z) / (p.x * d.x + p.y * d.y + p.z * d.z);
			// if behind the PhysicsMath::Triangle
			if (t >= 0.0f && t <= 1)
			{
				glm::vec3 intPoint = q1 + (d*t);
				// if the point of intersection is inside of the PhysicsMath::Triangle
				if (PhysicsMath::isInside(&intPoint, (&b2->tris[j])->t1, (&b2->tris[j])->t2, (&b2->tris[j])->t3))
				{
					t = PhysicsMath::pointPlaneDistance(&p, &pt);

					// update furthest point in
					if (t < furthest)
					{
						vert = &b1->verts[i];
						tri = &b2->tris[j];
						furthest = t;
						//VertPhysicsMath::Triangle vt;
						//vt.d = t;
						//vt.PhysicsMath::Triangle = &r.PhysicsMath::Triangles[j];
						//vt.vert = &verts[i];
						//vertPhysicsMath::Triangles.push_back(vt);
					}
				}
			}
		}
	}

	//for (int i = 0; i < vertTris.size(); ++i)
	{
		//vert = vertPhysicsMath::Triangles[i].vert;
		//PhysicsMath::Triangle = vertPhysicsMath::Triangles[i].PhysicsMath::Triangle;

		// collision response
		if (vert != nullptr && tri != nullptr)
		{
			//cout << "Collision!" << endl;

			// push from body r

			// normal of the PhysicsMath::Triangle that is bouncing off
			glm::vec3 n = glm::normalize(PhysicsMath::TriangleNormal(tri));
			// vector from furthest point to plane orthoganol
			glm::vec3 nt = n * furthest;

			if (!b1->isStatic && !b2->isStatic)
			{
				//cout << 1 << endl;
				b1->pos -= nt / 2.0f;
				b2->pos += nt / 2.0f;
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

			// collision point on PhysicsMath::Triangle/plane
			glm::vec3 p = (*vert + b1->pos);

			//cout << b1->vel.magnitude() << endl;

			// calc impulse
			float mu = 1;

			// collision point vectors
			glm::vec3 r1 = p - b1->pos;
			glm::vec3 r2 = p - b2->pos;

			// inertia tensors
			float * i1_1 = PhysicsMath::invert(glm::value_ptr(b1->I)); //r1 * b1->rotVec.w;//w1.normalized();
			float * i2_1 = PhysicsMath::invert(glm::value_ptr(b2->I));					   //float * i2_1 = invert(b2->I); // r2 * b2->rotVec.w;//w2.normalized();
			
			float * b1i = glm::value_ptr(b1->I);
			float i1 = glm::dot(glm::vec3(
				n.x * b1i[0] + n.y * b1i[3] + n.z * b1i[6],
				n.x * b1i[1] + n.y * b1i[4] + n.z * b1i[7],
				n.x * b1i[2] + n.y * b1i[5] + n.z * b1i[8]
			), n);
			float * b2i = glm::value_ptr(b2->I);
			float i2 = glm::dot(glm::vec3(
				n.x * b2i[0] + n.y * b2i[3] + n.z * b2i[6],
				n.x * b2i[1] + n.y * b2i[4] + n.z * b2i[7],
				n.x * b2i[2] + n.y * b2i[5] + n.z * b2i[8]
			), n);

			// relative vel of CoMs
			glm::vec3 vr = b2->vel - b1->vel;

			// impulse ratio

			//float jr =
			//	(-(1 + min(b1->e, b2->e)) * glm::dot(vr, n))
			//	/
			//	(1 / b1->mass + 1 / b2->mass +
			//	(glm::dot(glm::cross(r1, n), glm::cross(r1, n)) / i1) +
			//		(glm::dot(glm::cross(r2, n), glm::cross(r2, n)) / i2)); 



			float jr1 = (-(1 + min(b1->e, b2->e)) * glm::dot(vr, n));
			float jr2 = glm::dot(n, glm::cross(b1->I * glm::cross(r1, n), r1)) + glm::dot(n, glm::cross(b2->I * glm::cross(r2, n), r2));

			if (!b1->isStatic)
				jr2 += 1 / b1->mass;
			if (!b2->isStatic)
				jr2 += 1 / b2->mass;

			float jr = jr1 / jr2;
			jr *= mu;

			// impulse vectors
			glm::vec3 impulse1 = n * (-jr / b1->mass);
			glm::vec3 impulse2 = n * (jr / b2->mass);

			// apply impulse
			b1->vel += impulse1;
			b2->vel += impulse2;

			// angular vector of collision point
			glm::vec3 w1 = PhysicsMath::getRotatedVector(b1->rot, r1);
			//glm::vec4 w2 = getRotatedVector(b2->rotQuart, r2);

			// angular impulse vectors
			glm::vec3 W1 = w1 - (glm::cross(r1, n)* (jr / i1));

			//glm::vec4 W2 = w2 + (m33v(i2_1, r2.Cross(n)) * jr);

			// calc new quarternion for rotation and 
			// apply impulse magnitude to the angular speed
	/*		b1->rot.x = W1.x;
			b1->rot.y = W1.y;
			b1->rot.z = W1.z;
			b1->rot.b = W1.length() / r1.length();*/

			//W2 = glm::normalize(glm::cross(W2, r2));
			//b2->rotVec = ;
			//b2->rotVec.w = W2.magnitude() / r1.magnitude();

		}
	}
}