#ifndef TOOLS
#define TOOLS

//// Version 1.0

#include <GL\glut.h>
#include <math.h>
#include <algorithm> 
#include <locale.h> 
#include <string>
#include <sstream>
#include <thread>
#include <vector>
#include <iostream>
#include <exception>
#include <math.h>

using namespace std;

namespace Tools
{
	//// KEYSTROKE TOOLS - Hugh Nixon ////

	// KEY LIST //
	unsigned char			keyList[][3] =
	{
			{ 0x41, 'a', 'A' },
			{ 0x42, 'b', 'B' },
			{ 0x43, 'c', 'C' },
			{ 0x44, 'd', 'D' },
			{ 0x45, 'e', 'E' },
			{ 0x46, 'f', 'F' },
			{ 0x47, 'g', 'G' },
			{ 0x48, 'h', 'H' },
			{ 0x49, 'i', 'I' },
			{ 0x4A, 'j', 'J' },
			{ 0x4B, 'k', 'K' },
			{ 0x4C, 'l', 'L' },
			{ 0x4D, 'm', 'M' },
			{ 0x4E, 'n', 'N' },
			{ 0x4F, 'o', 'O' },
			{ 0x50, 'p', 'P' },
			{ 0x51, 'q', 'Q' },
			{ 0x52, 'r', 'R' },
			{ 0x53, 's', 'S' },
			{ 0x54, 't', 'T' },
			{ 0x55, 'u', 'U' },
			{ 0x56, 'v', 'V' },
			{ 0x57, 'w', 'W' },
			{ 0x58, 'x', 'X' },
			{ 0x59, 'y', 'Y' },
			{ 0x5A, 'z', 'Z' },
			{ 0x30, '0', ')' },
			{ 0x31, '1', '!' },
			{ 0x32, '2', '"' },
			{ 0x33, '3', '£' },
			{ 0x34, '4', '$' },
			{ 0x35, '5', '%' },
			{ 0x36, '6', '^' },
			{ 0x37, '7', '&' },
			{ 0x38, '8', '*' },
			{ 0x39, '9', '(' },
			{ 0xDC, '\\', '|' },
			{ 0xBC, ',', '<' },
			{ 0xBE, '.', '>' },
			{ 0xBF, '/', '?' },
			{ 0xBA, ';', ':' },
			{ 0xC0, '\'', '@' },
			{ 0xDE, '#', '~' },
			{ 0xBD, '-', '_' },
			{ 0xBB, '=', '+' },
			{ 0xDB, '[', '{' },
			{ 0xDD, ']', '}' },
			{ 0x20, ' ', ' ' }
	};

	// KEYBOARD STRUCT // 
	struct	CKeyboardController
	{
		char oldKeys[256];
		char newKeys[256];

		int focusOnMyWindow;
		string windowName;
		HWND win;

		CKeyboardController()
		{
			//LoadLibraryA("HUtils.dll");
			for (int i = 0; i < 256; i++)
			{
				oldKeys[i] = 0;
				newKeys[i] = 0;
			}
			focusOnMyWindow = 0;
			win = 0;
			windowName = "";
		}

		void run()
		{
			if (focusOnMyWindow)
			{
				if (GetForegroundWindow() != win)
				{
					for (int i = 0; i < 256; i++)
					{
						oldKeys[i] = 0;
						newKeys[i] = 0;
					}
					return;
				}
			}

			for (int i = 0; i < 256; i++)
			{
				oldKeys[i] = newKeys[i];
				newKeys[i] = !!GetAsyncKeyState(i);
			}
		}

		void setFocusWindow(string _windowName)
		{
			windowName = _windowName;
			win = FindWindowA(NULL, _windowName.c_str());
			focusOnMyWindow = 1;
		}

		vector<int> getAllDownEvent()
		{
			vector<int> v;
			for (int i = 0; i < 256; ++i)
				if (keyDownEvent(i))
					v.push_back(i);
			return v;
		}
		vector<int> getAllUpEvent()
		{
			vector<int> v;
			for (int i = 0; i < 256; ++i)
				if (keyUpEvent(i))
					v.push_back(i);
			return v;
		}
		vector<int> getAllDown()
		{
			vector<int> v;
			for (int i = 0; i < 256; ++i)
				if (keyDown(i))
					v.push_back(i);
			return v;
		}

		bool keyDown(int key)
		{
			return newKeys[key];
		}
		bool keyUpEvent(int key)
		{
			return !newKeys[key] && oldKeys[key];
		}
		bool keyDownEvent(int key)
		{
			return newKeys[key] && !oldKeys[key];
		}
		char getChar()
		{
			for (int i = 0; i < 48; ++i)
				if (keyDownEvent(keyList[i][0]))
					if (keyDown(VK_SHIFT))
						return keyList[i][2];
					else
						return keyList[i][1];
			return -1;
		}
	};
	CKeyboardController *	keyboard = new CKeyboardController();


	//// GL TOOLS ////

	static const float	pi = 3.14159265359f;
	static const float	pi2 = 6.28318530718f;

	HWND				win;

	// vector 4 struct //
	struct vec4
	{
		float x, y, z, w;

		vec4(float nx)
		{
			x = y = z = nx;
			w = 1;
		}
		vec4(float nx, float ny)
		{
			x = nx;
			y = ny;
			z = 0;
			w = 0;
		}
		vec4(float nx, float ny, float nz)
		{
			x = nx;
			y = ny;
			z = nz;
			w = 0;
		}
		vec4(float nx, float ny, float nz, float nw)
		{
			x = nx;
			y = ny;
			z = nz;
			w = nw;
		}
		vec4()
		{
			x = 0;
			y = 0;
			z = 0;
			w = 0;
		}

		vec4 &operator *	(vec4 a)
		{
			return vec4(a.x * x, a.y * y, a.z * z, w);
		}
		vec4 &operator *	(float a)
		{
			return vec4(a * x, a * y, a * z, w);
		}
		vec4 &operator *=	(vec4 a)
		{
			x *= a.x;
			y *= a.y;
			z *= a.z;
			return *this;
		}
		vec4 &operator *=	(float a)
		{
			x *= a;
			y *= a;
			z *= a;
			return *this;
		}

		vec4 &operator /	(float a)
		{
			return vec4(x / a, y / a, z / a, w);
		}
		vec4 &operator /=	(float a)
		{
			x /= a;
			y /= a;
			z /= a;
			return *this;
		}
		vec4 &operator /	(vec4 a)
		{
			return vec4(x / a.x, y / a.y, z / a.z, w);
		}
		vec4 &operator /=	(vec4 a)
		{
			x /= a.x;
			y /= a.y;
			z /= a.z;
			return *this;
		}

		vec4 &operator +	(vec4 a)
		{
			return vec4(a.x + x, a.y + y, a.z + z, w);
		}
		vec4 &operator +	(float a)
		{
			return vec4(a + x, a + y, a + z, w);
		}
		vec4 &operator +=	(vec4 a)
		{
			x += a.x;
			y += a.y;
			z += a.z;
			return *this;
		}
		vec4 &operator +=	(float a)
		{
			x += a;
			y += a;
			z += a;
			return *this;
		}

		vec4 &operator -	(vec4 a)
		{
			return vec4(x - a.x, y - a.y, z - a.z, w);
		}
		vec4 &operator -	(float a)
		{
			return vec4(x - a, y - a, z - a, w);
		}
		vec4 &operator -=	(vec4 a)
		{
			x -= a.x;
			y -= a.y;
			z -= a.z;
			return *this;
		}
		vec4 &operator -=	(float a)
		{
			x -= a;
			y -= a;
			z -= a;
			return *this;
		}
		
		bool operator ==	(vec4 a)
		{
			return x == a.x && y == a.y && z == a.z && w == a.w;
		}
		bool operator !=	(vec4 a)
		{
			return x != a.x && y != a.y && z != a.z && w != a.w;
		}

		inline vec4			inverse()
		{
			return vec4(1/x,1/y,1/z, w);
		}

		inline void			normalize()
		{
			float m = sqrtf(x * x + y * y + z * z);
			x /= m;
			y /= m;
			z /= m;
		}
		inline vec4			normalized()
		{
			float m = sqrtf(x * x + y * y + z * z);
			return vec4(x / m, y / m, z / m, w);
		}
		inline float		magnitude(vec4 in)
		{
			return sqrtf(in.x * in.x + in.y * in.y + in.z * in.z);
		}
		inline float		magnitude()
		{
			return sqrtf(x * x + y * y + z * z);
		}
		inline float		magnitude2()
		{
			return (x * x + y * y + z * z);
		}
		inline void			rotateZ(float r)
		{
			float currentRot = atan2(y, x);
			float mag = magnitude();
			x = mag * cos(currentRot + r);
			y = mag * sin(currentRot + r);
		}
		inline vec4			rotatedZ(float r)
		{
			float currentRot = atan2(y, x);
			float mag = magnitude();
			return vec4(mag * cos(currentRot + r), mag * sin(currentRot + r), z, w);
		}
		inline vec4			normalZ()
		{
			return rotatedZ(pi / 2);
		}
		inline void			rotateY(float r)
		{
			float currentRot = atan2(z, x);
			float mag = magnitude();
			x = mag * cos(currentRot + r);
			z = mag * sin(currentRot + r);
		}
		inline vec4			rotatedY(float r)
		{
			float currentRot = atan2(z, x);
			float mag = magnitude();
			return vec4(mag * cos(currentRot + r), y, mag * sin(currentRot + r), w);
		}
		inline vec4			normalY()
		{
			return rotatedY(pi / 2);
		}
		inline void			rotate(float angler, vec4 k)
		{
			// Rodrigues formula
			k.normalize();
			vec4 v = *this;
			*this += v * cos(angler) + k * (1 - cos(angler)) * k.Dot(v) + k.Cross(v) * sin(angler);
		}
		inline vec4			rotated(float angler, vec4 k)
		{
			k.normalize();
			vec4 v = *this;
			return *this + v * cos(angler) + k * (1 - cos(angler)) * k.Dot(v) + k.Cross(v) * sin(angler);
		}
		inline float		dot(vec4 a)
		{
			return x * a.x + y * a.y + z * a.z;
		}
		inline float		Dot(vec4 a)
		{
			return (x * a.x + y * a.y + z * a.z) / (magnitude() * a.magnitude());
		}
		inline vec4			Cross(vec4 a)
		{
			return vec4(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
		}
		inline void			glColour4()
		{
			glColor4f(x, y, z, w);
		}
		inline void			glColour3()
		{
			glColor3f(x, y, z);
		}
		inline std::string	toString()
		{
			std::stringstream s;
			s << "| X: " << x << " | Y: " << y << " | Z: " << z << " | W: " << w;
			return s.str();
		}
		inline void			print()
		{
			cout << toString() << endl;
		}
		inline void			glVertex3()
		{
			glVertex3f(x, y, z);
		}
		inline void			glVertex2()
		{
			glVertex2f(x, y);
		}
		inline void			glTranslate()
		{
			glTranslatef(x, y, z);
		}
		inline vec4			setAlpha(float a)
		{
			return vec4(x, y, z, a);
		}

		inline float		qMag()
		{
			return sqrtf(x * x + y * y + z * z + w * w);
		}
		inline vec4			qNormalized()
		{
			float d = qMag();
			return vec4(x / d, y / d, z / d, w / d);
		}
		inline vec4			qNormalize()
		{
			float d = qMag();
			x /= d;
			y /= d;
			z /= d;
			w /= d;
		}
		inline vec4			getConQuart()
		{
			vec4 t = normalized();
			float alpha = sinf(t.w / 2.0f);
			return vec4(t.x * alpha, t.y * alpha, t.z * alpha, cosf(t.w / 2.0f));
		}
		inline vec4			conQuart()
		{
			normalize();
			float alpha = sinf(w / 2.0f);
			x = x * alpha;
			y = y * alpha;
			z = z * alpha;
			w = cosf(w / 2.0f);
			return *this;
		}
		inline vec4			getQConj()
		{
			return vec4(-x, -y, -z, w);
		}
		inline vec4			qMult(vec4 r)
		{
			vec4 t;
			t.w = r.w*w - r.x*x - r.y*y - r.z*z;
			t.x = r.w*x + r.x*w - r.y*z + r.z*y;
			t.y = r.w*y + r.x*z + r.y*w - r.z*x;
			t.z = r.w*z - r.x*y + r.y*x + r.z*w;
			return t;
		}
		inline float*		getRotationMatrix()
		{
			vec4 q = normalized();
			float * mat = new float[9];
			mat[0] = 1 - 2 * (q.y*q.y + q.z*q.z);
			mat[1] = 2 * (q.x*q.y - q.z*q.w);
			mat[2] = 2 * (q.x*q.z + q.y*q.w);
			mat[3] = 2  * (q.x*q.y + q.z*q.w);
			mat[4] = 1 - 2 * (q.x*q.x + q.z*q.z);
			mat[5] = 2 * (q.y*q.z + q.x*q.z);
			mat[6] = 2 * (q.x*q.z - q.y*q.z);
			mat[7] = 2  * (q.y*q.z + q.x*q.z);
			mat[8] = 1 - 2 * (q.x*q.x + q.y*q.y);
			return mat;
		}
		inline float*		getRotationMatrixOrtho()
		{
			vec4 q = normalized();
			float * mat = new float[9];
			mat[0] = w*w + x*x - y*y - z*z;
			mat[1] = 2 * x*y - 2 * w*z;
			mat[2] = 2 * x*z + 2 * w*y;
			mat[3] = 2 * x*y + 2 * w*z;
			mat[4] = w*w - x*x + y*y - z*z;
			mat[5] = 2 * y*z - 2 * w*x;
			mat[6] = 2 * x*z - 2 * w*y;
			mat[7] = 2 * y*z + 2 * w*x;
			mat[8] = w*w - x*x - y*y + z*z;
			return mat;
		}
		inline vec4			getRotatedPoint(vec4 p)
		{
			vec4 r = HP(HP(*this, p), (*this).getQConj());
			r.w = 0;
			return r;
		}
		inline void			rotateArray(vec4 * arr, int n)
		{
			for (int i = 0; i < n; ++i)
				arr[i] = getRotatedPoint(arr[i]);
		}
		inline void			rotateArray(vector<vec4> * arr)
		{
			for (int i = 0; i < arr->size(); ++i)
				(*arr)[i] = getRotatedPoint((*arr)[i]);
		}
		// Hamilton Procduct of two quarternions
		static vec4			HP(vec4 _1, vec4 _2)
		{
			vec4 r;
			r.w = _1.w*_2.w - _1.x*_2.x - _1.y*_2.y - _1.z*_2.z;
			r.x = _1.w*_2.x + _1.x*_2.w + _1.y*_2.z - _1.z*_2.y;
			r.y = _1.w*_2.y - _1.x*_2.z + _1.y*_2.w + _1.z*_2.x;
			r.z = _1.w*_2.z + _1.x*_2.y - _1.y*_2.x + _1.z*_2.w;
			return r;
		}
	};

	// MATRIX 4 //
	struct mat4
	{
		float m[16];

		mat4()
		{
			for (int i = 0; i < 16; ++i)
				m[i] = 0;
		}
		mat4(float * m)
		{
			for (int i = 0; i < 16; ++i)
				this->m[i] = m[i];
		}
	};

	// TRI //
	struct tri
	{
		vec4 *t1, *t2, *t3;
		tri() {}
		tri(vec4* t1, vec4* t2, vec4* t3)
		{
			this->t1 = t1;
			this->t2 = t2;
			this->t3 = t3;
		}
		void draw()
		{
			t1->glVertex3();
			t2->glVertex3();
			t3->glVertex3();
		}
	};

	// timer struct // 
	struct Timer
	{
		float init = 0, current = 0;
		bool cyclic = 0;
		Timer()
		{}
		Timer(float i)
		{
			init = current = i;
		}
		Timer(float i, bool c)
		{
			init = current = i;
			cyclic = c;
		}
		Timer(float i, bool c, bool z)
		{
			init = current = i;
			cyclic = c;
			if (z)
				current = 1;
		}
		inline bool tick()
		{
			if (current > 0)
				current--;
			else if (cyclic)
				reset();
			return zero();
		}
		inline void reset()
		{
			current = init;
		}
		inline bool zero()
		{
			return current == 0;
		}
	};


	// MISC FUNCTIONS //
	inline float		roundf(float f, int c)
	{

		return ((float)((int)(f * c)) / c);
	}
	inline char *		commaprint(unsigned long n)
	{
		static int comma = '\0';
		static char retbuf[30];
		char *p = &retbuf[sizeof(retbuf) - 1];
		int i = 0;

		if (comma == '\0') {
			struct lconv *lcp = localeconv();
			if (lcp != NULL) {
				if (lcp->thousands_sep != NULL &&
					*lcp->thousands_sep != '\0')
					comma = *lcp->thousands_sep;
				else	comma = ',';
			}
		}

		*p = '\0';

		do {
			if (i % 3 == 0 && i != 0)
				*--p = comma;
			*--p = '0' + n % 10;
			n /= 10;
			i++;
		} while (n != 0);

		return p;
	}
	inline float		randf()
	{
		return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	}
	inline vec4			randv()
	{
		return vec4(randf(), randf(), randf(), randf());
	}
	inline float		radDeg(float r)
	{
		return r / pi * 180;
	}
	inline float		degRad(float d)
	{
		return d / 180 * pi;
	}
	inline std::string	repeatChar(char c, int amount)
	{
		std::string str = "";
		for (int i = 0; i < amount; ++i)
			str += c;
		return str;
	}
	inline bool			isNumber(const std::string& s)
	{
		return !s.empty() && std::find_if(s.begin(),
			s.end(), [](char c) { return !isdigit(c); }) == s.end();
	}

	// MATH FUNCTIONS //
	inline bool			intersect(vec4 a0, vec4 a1, vec4 b0, vec4 b1)
	{
		vec4 a = (a1 - a0) * -1;
		vec4 b = (b1 - b0) * -1;
		float mu = ((b0.y - a0.y)*a.x / a.y + a0.x - b0.x) / (b.x - b.y*a.x / a.y);
		float lambda = ((b0.x - a0.x) + mu*b.x) / a.x;
		bool p1 = (abs((a0.x + lambda * a.x) - (b0.x + mu * b.x)) < 0.001f && abs((a0.y + lambda * a.y) - (b0.y + mu * b.y)) < 0.01f);
		bool p2 = lambda > 0.0f && lambda < 1.0f;
		bool p3 = mu > 0.0f && mu < 1.0f;
		return p1 && p2 && p3;
	}
	inline float		lineDistance(vec4 p0, vec4 p1, vec4 p2)
	{
		vec4 p0p1 = p0 - p1;
		vec4 p1p2 = p2 - p1;
		return (p1p2.x*p0p1.y - p0p1.x*p1p2.y) / (p1p2.magnitude());
	}
	inline vec4			getProjectedPointLine(vec4 p0, vec4 p1, vec4 p2)
	{
		vec4 p0p1 = p0 - p1;
		vec4 p1p2 = p2 - p1;
		float dp = p0p1.Dot(p1p2);
		float len1 = p1p2.magnitude();
		float len2 = p0p1.magnitude();
		float c = dp / (len1*len2);
		float t = c * len2;
		return p1 + (p1p2 * (t / len1));
	}
	inline float		getProjectedT(vec4 p0, vec4 p1, vec4 p2)
	{
		vec4 p0p1 = p0 - p1;
		vec4 p1p2 = p2 - p1;
		float dp = p0p1.Dot(p1p2);
		float len1 = p1p2.magnitude();
		float len2 = p0p1.magnitude();
		float c = dp / (len1*len2);
		float t = c * len2;
		return t / len1;
	}
	inline bool			isPointWithinLine(float distance, vec4 p0, vec4 p1, vec4 p2)
	{
		float t = getProjectedT(p0, p1, p2);
		return abs(lineDistance(p0, p1, p2)) < distance && t > 0.0f && t < 1.0f;
	}
	inline float		triArea3(vec4 p0, vec4 p1, vec4 p2)
	{
		float d1 = (p1.y - p0.y)*(p2.z - p0.z) - (p1.z - p0.z)*(p2.y - p0.y);
		float d2 = (p1.z - p0.z)*(p2.x - p0.x) - (p1.x - p0.x)*(p2.z - p0.z);
		float d3 = (p1.x - p0.x)*(p2.y - p0.y) - (p1.y - p0.y)*(p2.x - p0.x);
		return sqrt(d1*d1 + d2*d2 + d3*d3) / 2;
	}
	inline bool			isInside(vec4 p, vec4 t1, vec4 t2, vec4 t3)
	{
		float totArea = triArea3(p, t1, t2) + triArea3(p, t2, t3) + triArea3(p, t3, t1);
		float actArea = triArea3(t1, t2, t3);
		return totArea <= actArea*1.005f;
	}
	inline float		triArea3(vec4 *p0, vec4 *p1, vec4 *p2)
	{
		float d1 = (p1->y - p0->y)*(p2->z - p0->z) - (p1->z - p0->z)*(p2->y - p0->y);
		float d2 = (p1->z - p0->z)*(p2->x - p0->x) - (p1->x - p0->x)*(p2->z - p0->z);
		float d3 = (p1->x - p0->x)*(p2->y - p0->y) - (p1->y - p0->y)*(p2->x - p0->x);

		return sqrtf(d1*d1 + d2*d2 + d3*d3) / 2;
	}
	inline bool			isInside(vec4 *p, vec4 *t1, vec4 *t2, vec4*t3)
	{
		float totArea = triArea3(p, t1, t2) + triArea3(p, t2, t3) + triArea3(p, t3, t1);
		float actArea = triArea3(t1, t2, t3);
		return totArea <= actArea*1.005f;
	}
	inline bool			isInsideRec(vec4 pos, vec4 size, vec4 p)
	{
		return 
			p.x >= pos.x && p.x <= pos.x + size.x && 
			p.y >= pos.y && p.y <= pos.y + size.y;
	}


	// GL DRAW FUNCTIONS //
	inline void			drawArrow()
	{
		glBegin(GL_LINE_STRIP);
		glVertex2f(1, 0);
		glVertex2f(0, 0);

		glVertex2f(0, 0);
		glVertex2f(0.2f, 0.1f);

		glVertex2f(0.2f, 0.1f);
		glVertex2f(0.2f, -0.1f);

		glVertex2f(0.2f, -0.1f);
		glVertex2f(0, 0);
		glEnd();
	}
	inline void			drawBullet()
	{
		glBegin(GL_LINE_STRIP);
		glVertex2f(0, 0.2f);
		glVertex2f(1, 0);

		glVertex2f(1, 0);
		glVertex2f(0, -0.2f);
		glEnd();
	}
	inline void			lineCirc(int iter)
	{
		// uses polar coords
		float ang = pi2 / float(iter);
		glBegin(GL_LINE_LOOP);
		for (int i = 0; i < iter; ++i)
			glVertex2f(cos(i*ang), sin(i*ang));
		glEnd();

	}
	inline void			fillCirc(int iter)
	{
		// uses polar coords
		float ang = pi2 / float(iter);
		glBegin(GL_POLYGON);
		for (int i = 0; i < iter; ++i)
			glVertex2f(cos(i*ang), sin(i*ang));
		glEnd();

	}
	inline void			fillRect()
	{
		// unit square matrix
		glBegin(GL_POLYGON);
		glVertex2f(0.0f, 0.0f);
		glVertex2f(0.0f, 1.0f);
		glVertex2f(1.0f, 1.0f);
		glVertex2f(1.0f, 0.0f);
		glEnd();
	}
	inline void			lineRect()
	{
		// unit square matrix
		glBegin(GL_LINE_LOOP);
		glVertex2f(0.0f, 0.0f);
		glVertex2f(0.0f, 1.0f);
		glVertex2f(1.0f, 1.0f);
		glVertex2f(1.0f, 0.0f);
		glEnd();
	}
	inline void			fillRRect(float f, int iter)
	{
		// draws four circles in the corners
		float ff = 1.0f - f;
		glPushMatrix();
		glTranslatef(0.5f, 0.5f, 0.0f);
		glScalef(0.5f, 0.5f, 0.0f);
		glPushMatrix();
		glTranslatef(-f, f, 0.0f);
		glScalef(ff, ff, 0.0f);
		fillCirc(iter);
		glPopMatrix();

		glPushMatrix();
		glTranslatef(f, f, 0.0f);
		glScalef(ff, ff, 0.0f);
		fillCirc(iter);
		glPopMatrix();

		glPushMatrix();
		glTranslatef(f, -f, 0.0f);
		glScalef(ff, ff, 0.0f);
		fillCirc(iter);
		glPopMatrix();

		glPushMatrix();
		glTranslatef(-f, -f, 0.0f);
		glScalef(ff, ff, 0.0f);
		fillCirc(iter);
		glPopMatrix();

		// draw overlapping rectangles
		//horr
		glPushMatrix();
		glScalef(2.0f, f*2.0f, 0.0f);
		glTranslatef(-0.5f, -0.5f, 0.0f);
		fillRect();
		glPopMatrix();
		//vert
		glPushMatrix();
		glScalef(f*2.0f, 2.0f, 0.0f);
		glTranslatef(-0.5f, -0.5f, 0.0f);
		fillRect();
		glPopMatrix();
		glPopMatrix();
	}
	inline void			vGRect(vec4 c1, vec4 c2)
	{
		// unit square with colours
		glBegin(GL_POLYGON);
		glColor4f(c1.x, c1.y, c1.z, c1.w);
		glVertex2f(0.0f, 0.0f);
		glColor4f(c2.x, c2.y, c2.z, c2.w);
		glVertex2f(0.0f, 1.0f);
		glColor4f(c2.x, c2.y, c2.z, c2.w);
		glVertex2f(1.0f, 1.0f);
		glColor4f(c1.x, c1.y, c1.z, c1.w);
		glVertex2f(1.0f, 0.0f);
		glEnd();
	}
	inline void			hGRect(vec4 c1, vec4 c2)
	{
		// unit square with colours
		glBegin(GL_POLYGON);
		glColor4f(c1.x, c1.y, c1.z, c1.w);
		glVertex2f(0.0f, 0.0f);
		glColor4f(c1.x, c1.y, c1.z, c1.w);
		glVertex2f(0.0f, 1.0f);
		glColor4f(c2.x, c2.y, c2.z, c2.w);
		glVertex2f(1.0f, 1.0f);
		glColor4f(c2.x, c2.y, c2.z, c2.w);
		glVertex2f(1.0f, 0.0f);
		glEnd();
	}
	inline void			drawText10(float x, float y, char ch)
	{
		// draw single character
		glRasterPos2f(x + 2.0f, y + 8.0f);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, ch);
	}
	inline void			drawText10(float x, float y, std::string text)
	{
		// draw string
		glRasterPos2f(x + 2.0f, y + 8.0f);
		for (int i = 0; i < text.length(); i++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, text[i]);
	}
	inline void			drawText12(float x, float y, std::string text)
	{
		glRasterPos2f(x + 2.0f, y + 16.0f);
		for (int i = 0; i < text.length(); i++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, text[i]);
	}
	inline void			drawText18(float x, float y, std::string text)
	{
		glRasterPos2f(x + 2.0f, y + 16.0f);
		for (int i = 0; i < text.length(); i++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
	}
	inline void			drawText24(float x, float y, std::string text)
	{
		glRasterPos2f(x + 2.0f, y + 16.0f);
		for (int i = 0; i < text.length(); i++)
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, text[i]);
	}
	inline void			drawTextX(float x, float y, std::string text, void * FONT)
	{
		glRasterPos2f(x, y);
		for (int i = 0; i < text.length(); i++)
			glutBitmapCharacter(FONT, text[i]);
	}
	inline void			drawTextX(float x, float y, char ch, void * FONT)
	{
		// draw single character
		glRasterPos2f(x , y);
		glutBitmapCharacter(FONT, ch);
	}
	inline void			lineCube()
	{
		glBegin(GL_LINE_LOOP);

		glVertex3f(0, 0, 0);

		glVertex3f(1, 0, 0);
		glVertex3f(1, 1, 0);
		glVertex3f(0, 1, 0);
		glVertex3f(0, 0, 0);

		glVertex3f(0, 0, 1);

		glVertex3f(1, 0, 1);
		glVertex3f(1, 0, 0);
		glVertex3f(1, 0, 1);

		glVertex3f(1, 1, 1);
		glVertex3f(1, 1, 0);
		glVertex3f(1, 1, 1);

		glVertex3f(0, 1, 1);
		glVertex3f(0, 1, 0);
		glVertex3f(0, 1, 1);

		glVertex3f(0, 0, 1);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1);

		glEnd();
	}
	inline void			lineSphere(int lats, int longs) 
	{
		int i, j;
		for (i = 0; i <= lats; i++)
		{
			float lat0 = pi * (-0.5 + (float)(i - 1) / lats);
			float z0 = sin(lat0);
			float zr0 = cos(lat0);

			float lat1 = pi * (-0.5 + (float)i / lats);
			float z1 = sin(lat1);
			float zr1 = cos(lat1);

			glBegin(GL_LINE_LOOP);
			for (j = 0; j <= longs; j++)
			{
				float lng = pi2 * (float)(j - 1) / longs;
				float x = cos(lng);
				float y = sin(lng);

				glNormal3f(x * zr0, y * zr0, z0);
				glVertex3f(x * zr0, y * zr0, z0);
				glNormal3f(x * zr1, y * zr1, z1);
				glVertex3f(x * zr1, y * zr1, z1);

			}
			glEnd();
		}
	}
	inline void			drawRecTex(float x, float y, float w, float h, GLuint tex, bool ver, bool hor)
	{
		glEnable(GL_TEXTURE_2D);

		glBindTexture(GL_TEXTURE_2D, tex);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glBegin(GL_QUADS);
		glColor3f(1, 1, 1);
		glTexCoord2f(0, 0);
		glVertex2f(x, y);
		glTexCoord2f(1 - 2 * hor, 0);
		glVertex2f(x + w, y);
		glTexCoord2f(1 - 2 * hor, 1 - 2 * ver);
		glVertex2f(x + w, y + h);
		glTexCoord2f(0, 1 - 2 * ver);
		glVertex2f(x, y + h);
		glEnd();

		glDisable(GL_TEXTURE_2D);
	}
	inline void			drawRecTex(float x, float y, float w, float h, GLuint tex)
	{
		glEnable(GL_TEXTURE_2D);

		glBindTexture(GL_TEXTURE_2D, tex);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		glBegin(GL_QUADS);
		glColor3f(1, 1, 1);
		glTexCoord2f(0, 0);
		glVertex2f(x, y);
		glTexCoord2f(1, 0);
		glVertex2f(x + w, y);
		glTexCoord2f(1, 1);
		glVertex2f(x + w, y + h);
		glTexCoord2f(0, 1);
		glVertex2f(x, y + h);
		glEnd();

		glDisable(GL_TEXTURE_2D);
	}

	// GL SHORTER FUNCTION NAMES //
	inline void push()
	{
		glPushMatrix();
	}
	inline void pop()
	{
		glPopMatrix();
	}
	inline void poly()
	{
		glBegin(GL_POLYGON);
	}
	inline void ge()
	{
		glEnd();
	}

	// GL MATRICIES //
	GLdouble		projection[16];
	GLdouble		modelview[16];
	GLint			viewport[4];
	inline void		updateGLMatrices()
	{
		glGetIntegerv(GL_VIEWPORT, viewport);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	}
	inline vec4		w2s(vec4 p)
	{
		GLdouble oX; GLdouble oY; GLdouble oZ;
		gluProject(p.x, p.y, p.z, modelview, projection, viewport, &oX, &oY, &oZ);
		return vec4(oX, viewport[3] - oY, 0);
	}

	// FPS COUNTERS //
	int				tfc = 0;
	int				ffc = 0;
	float			tps = 0;
	float			fps = 0;
	inline void		incFFC()
	{
		ffc++;
	}
	void			fpsThread()
	{
		while (1)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(500));
			fps = ffc * 2;
			tps = tfc * 2;
			ffc = tfc = 0;
		}
	}
	void			detachFPSThread()
	{
		std::thread fps(fpsThread);
		fps.detach();
	}
	void			drawFPS()
	{
		glPushMatrix();
		glTranslatef(0, 18, 0);
		std::stringstream s;
		s << "FPS " << fps;
		drawText10(0, 0, s.str());
		glPopMatrix();

		glPushMatrix();
		glTranslatef(0, 18 * 2, 0);
		s.str(std::string());
		s << "TPS " << tps;
		drawText10(0, 0, s.str());
		glPopMatrix();
	}

	// OTHER GL STUFF //
	vec4	dir2mouse		= vec4();
	vec4	mousePos		= vec4();
	vec4	mouseDir		= vec4();
	vec4	window_i		(400, 400);
	vec4	window			= window_i;
	vec4	windowPos		(0, 0);
	string	windowName		= "lel";
	float	fov				= 60;
	vec4	worldPosition	(0, 0, -5);
	vec4	lookDirection	(0, 0, 1);
	vec4	upVec			(0, 1, 0);
	vec4	clearColour		(0, 0, 0);
	float	mouseSensitivity = 0.01f;
	bool	mouseLock		= false;
	int		constWaitTime	= 10;

	inline void start3dDraw()
	{
		glEnable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		vec4 lookPosition = lookDirection + worldPosition;
		gluLookAt(worldPosition.x, worldPosition.y, worldPosition.z, lookPosition.x, lookPosition.y, lookPosition.z, 0, 1, 0);
		updateGLMatrices();
	}
	inline void start2dDraw()
	{
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glOrtho(0.0, window.x, window.y, 0.0, -1.0, 10.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
	inline void end2dDraw()
	{
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
	}

	std::chrono::high_resolution_clock::time_point start, finish;
	inline void startProcessLoop()
	{
		start = std::chrono::high_resolution_clock::now();
		POINT cursorPos;
		GetCursorPos(&cursorPos);
		vec4 oldMouse = mousePos;
		mousePos = vec4(cursorPos.x - glutGet(GLUT_WINDOW_X), cursorPos.y - glutGet(GLUT_WINDOW_Y));
		mouseDir = oldMouse - mousePos;
		dir2mouse = vec4(window.x / 2, window.y / 2) - mousePos;

		if (mouseLock)
		{
			glutWarpPointer(window.x / 2, window.y / 2);
			glutSetCursor(GLUT_CURSOR_NONE);
		}
		else
			glutSetCursor(GLUT_CURSOR_INHERIT);
	}
	inline void endProcessLoop()
	{
		finish = std::chrono::high_resolution_clock::now();
		int ms = float(std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count());
		long newWait = constWaitTime - ms;
		newWait = newWait < 0 ? 0 : newWait;
		std::this_thread::sleep_for(std::chrono::milliseconds(newWait));
		tfc++;
	}

	inline void fpsCamera()
	{
		float pitch = dir2mouse.y * mouseSensitivity;
		float yaw = dir2mouse.x * mouseSensitivity;
		vec4 q = upVec.Cross(lookDirection);
		lookDirection.rotateY(yaw);
		lookDirection.rotate(pitch, q);		
		lookDirection.normalize();		
	}

	//// GUI TOOLS ////

	// constants //
	static HCURSOR		hCursorArrow	= LoadCursor(NULL, IDC_ARROW);
	static HCURSOR		hCursorText		= LoadCursor(NULL, IDC_IBEAM);

	// component bitmask flags //
	const unsigned int	MOUSE_DOWN						= 1;
	const unsigned int	COMPONENT_TOGGLABLE				= 2;
	const unsigned int	COMPONENT_TOGGLED				= 4;
	const unsigned int	COMPONENT_VISIBLE				= 8;
	const unsigned int	COMPONENT_ENABLED				= 16;
	const unsigned int	COMPONENT_PRESSED				= 32;
	const unsigned int	COMPONENT_FOCUSED				= 64;
	const unsigned int	COMPONENT_MULTILINE				= 128;
	const unsigned int	COMPONENT_BORDER				= 256;
	const unsigned int	COMPONENT_ROUNDED				= 512;
	const unsigned int	COMPONENT_VERTICAL_GRADIENT		= 1024;
	const unsigned int	COMPONENT_HORIZONTAL_GRADIENT	= 2048;
	const unsigned int	COMPONENT_REAL_DELTA			= 4096;
	const unsigned int	COMPONENT_CENSORED				= 8192;
	const unsigned int	COMPONENT_OPAQUE				= 16384;
	const unsigned int	COMPONENT_CENTERED				= 32768;
	const unsigned int	COMPONENT_LEFT					= 65536;

	// predefined colours // 
	vec4 COLOUR_WHITE		(1.0f, 1.0f, 1.0f, 1);
	vec4 COLOUR_BLACK		(0.0f, 0.0f, 0.0f, 1);
	vec4 COLOUR_LIGHT_GREY	(0.8f, 0.8f, 0.8f, 1);
	vec4 COLOUR_DARK_GREY	(0.3f, 0.3f, 0.3f, 1);
	vec4 COLOUR_GREY		(0.5f, 0.5f, 0.5f, 1);
	vec4 COLOUR_CYAN		(0.0f, 1.0f, 1.0f, 1);
	vec4 COLOUR_RED			(1.0f, 0.0f, 0.0f, 1);
	vec4 COLOUR_GREEN		(0.0f, 1.0f, 0.0f, 1);
	vec4 COLOUR_PURPLE		(1, 0, 1, 1);
	vec4 COLOUR_BLUE		(0.0f, 0.0f, 1.0f, 1);
	vec4 COLOUR_YELLOW		(1, 1, 0, 1);
	vec4 COLOUR_OFF_WHITE	(0.96f, 0.96f, 0.96f, 1);
	vec4 COLOUR_OFF_BLACK	(0.18f, 0.18f, 0.18f, 1);
	vec4 COLOUR_ORANGE		(1, 0.5, 0, 1);

	// prototypes //
	class GUIManager;
	class Container;
	class Component;
	class Button;
	class TextField;
	class Label;
	class Slider;

	// classes //
	class Component
	{
	public:
		Component()
		{
			text = downText = "";
			downEvent = upEvent = moveEvent = nullptr;
			flags = COMPONENT_ENABLED | COMPONENT_VISIBLE;
		}

		GUIManager * gui;
		Container * container;
		vec4 pos, size;
		vec4 bColour, fColour, textColour;
		string text, downText;
		void * FONT = GLUT_BITMAP_8_BY_13;
		void(*downEvent)();
		void(*upEvent)();
		void(*moveEvent)();
		int id = -1;
		unsigned int flags;
		int componentUID = 0;

		virtual void draw() = 0;
		virtual void handleEvents() = 0;
		bool isInside();
		bool test(unsigned int f);
		void setText(string text);
		void appendText(string text);
	};
	class Container : public Component
	{
	public:
		vector < Component * > components;

		void add(Component *c);
		void draw();
		void handleEvents();
		Component * getComponent(int _id);

		Container() : Component()
		{
			componentUID = 1;
		}
		Container(unsigned int f) : Component()
		{
			flags |= f;
		}

	};
	class Button : public Component
	{
	public:
		void draw();
		void handleEvents();

		Button() : Component()
		{
			componentUID = 2;
		}
		Button(void dE(), void uE(), unsigned int f) : Button()
		{
			downEvent = dE; upEvent = uE;
			flags |= f;
		}

		void setText(string txt);
		void setText(string upTxt, string downTxt);

	private:
	};
	class TextField : public Component
	{
	public:
		int cursorPos, startPos, endPos;
		TextField * tabTo = nullptr;
		TextField * tabBack = nullptr;

		void handleEvents();

		void draw();

		int findCursorPos();

		void getText(string * str)
		{
			str = &text;
		}

		TextField() : Component()
		{
			componentUID = 3;
			cursorPos = 0;
			startPos = 0;
			endPos = 0;
		}

		TextField(unsigned int f) : TextField()
		{
			flags |= f;
		}
	};
	class Label : public Component
	{
	public:
		int startCursor = 0;
		int endCursor = 0;
		int lines = 0;

		void draw();
		void handleEvents();
		void setTextMode(int n);
		void setText(string _text);
		void appendText(string text);
		void scrollUp();
		void scrollDown();

		Label() : Component()
		{
			componentUID = 4;
			textMode = 0;
		}
		Label(unsigned int f) : Label()
		{

			flags |= f;
			textMode = 0;
		}

	private:
		int textMode;
	};
	class Slider : public Component
	{
	public:

		int min, max, value;
		bool vertical = false;
		bool alreadyDown = false;
		int direction = 0;

		void draw();
		void handleEvents();
		bool isInsideEx();
		void setValue(int v);

		Slider() : Component()
		{
			min = 0; max = 10; value = 0;
			componentUID = 5;
		}
		Slider(void mE(), unsigned int f) : Slider()
		{
			moveEvent = mE;
			flags |= f;
		}
	};
	class GUIManager
	{
	public:
		bool locked = false;
		int id_held = -1;
		int ids = 0;
		vector<char> charBuffer;
		vector < Component * > components;
		float blink = 0.0f;;
		POINT cursorPos;
		vec4 mousePos, mouseOffset;
		bool moving = false;
		bool moved = false;
		Component * from = nullptr;
		Component * to = nullptr;
		int cDown = 0;
		volatile bool mouseOnGUI = 0;
		volatile int amountOn = 0;
		string windowName = "Default Name";
		Container * background = new Container();

		GUIManager()
		{
			background->pos = vec4();
			background->size = vec4(4096, 4096);
			background->container = nullptr;
		}

		void add(Component *c);
		void draw();
		void handleEvents();
		void updateMousePosition();
		void windowMove();
		void focusChange();
		void updateFocus(Component * _from, Component * _to);
		void changeSize(float w, float h);
		void init();
		void removeWindowBorder();
		void minimise();
		void terminate();
		void setDimensions(float w, float h);
	private:

	}  *gui = new GUIManager();

	// component - super class
	bool Component::isInside()
	{
		POINT mp;
		GetCursorPos(&mp);
		mp.x -= glutGet(GLUT_WINDOW_X);
		//mp.x *= gui->window_w_i / gui->window_w;
		mp.y -= glutGet(GLUT_WINDOW_Y);
		//mp.y *= gui->window_h_i / gui->window_h;
		if (container != nullptr)
		{
			if (mp.x >= pos.x + container->pos.x + (test(COMPONENT_CENTERED)*(+container->size.x / 2 - size.x / 2))
				&& mp.x <= pos.x + size.x + container->pos.x + (test(COMPONENT_CENTERED)*(+container->size.x / 2 - size.x / 2))
				&& mp.y >= pos.y + container->pos.y
				&& mp.y <= pos.y + size.y + container->pos.y)
				return true;
			else
				return false;
		}
		else
			if (mp.x >= pos.x
				&& mp.x <= pos.x + size.x
				&& mp.y >= pos.y
				&& mp.y <= pos.y + size.y)
				return true;
			else
				return false;
	}
	bool Component::test(unsigned int f)
	{
		return (flags & f) == f;
	}
	void Component::setText(string text)
	{
	}

	// container
	void Container::handleEvents()
	{
		if (test(COMPONENT_VISIBLE) && test(COMPONENT_ENABLED))
			for (int i = 0; i < components.size(); ++i)
				components[i]->handleEvents();
	}
	void Container::draw()
	{
		if (test(COMPONENT_VISIBLE))
		{
			bool opaque = test(COMPONENT_OPAQUE);
			// border
			if (test(COMPONENT_BORDER))
			{
				glPushMatrix();
				glColor4f(fColour.x, fColour.y, fColour.z, fColour.w + opaque);
				glTranslatef(pos.x - 1, pos.y - 1, 0.0f);
				glScalef(size.x + 2, size.y + 2, 0.0f);
				if (test(COMPONENT_ROUNDED))
					fillRRect(0.9f, 10);
				else
					fillRect();
				glPopMatrix();
			}

			// inside
			glPushMatrix();
			glColor4f(bColour.x, bColour.y, bColour.z, bColour.w + opaque);
			glTranslatef(pos.x, pos.y, 0.0f);
			glScalef(size.x, size.y, 0.0f);
			if (test(COMPONENT_ROUNDED))
				fillRRect(0.9f, 10);
			else
			{
				if (test(COMPONENT_VERTICAL_GRADIENT))
					vGRect(fColour, bColour);
				else if (test(COMPONENT_HORIZONTAL_GRADIENT))
					hGRect(fColour, bColour);
				else
					fillRect();
			}
			glPopMatrix();

			// draw components it holds
			for (int i = components.size() - 1; i >= 0; --i)
				components[i]->draw();
		}

	}
	void Container::add(Component *c)
	{
		c->container = this;
		c->id = gui->ids++;
		c->gui = gui;
		components.push_back(c);
	}

	// slider
	void Slider::setValue(int v)
	{
		// clamp
		value = v < 0 ? 0 : v > max ? max : v;
	}
	void Slider::draw()
	{
		if (container != nullptr)
		{
			if (test(COMPONENT_VISIBLE))
			{
				bool opaque = test(COMPONENT_OPAQUE);
				bool isDown = test(COMPONENT_PRESSED);
				vec4 c1 = bColour;
				vec4 c2 = fColour;
				// text inbevel amount
				float shift = 0.0f;
				// invert colours if down
				if (isDown)
				{
					c1 = vec4(1.0f - c1.x, 1.0f - c1.y, 1.0f - c1.z, c1.w + opaque);
					c2 = vec4(1.0f - c2.x, 1.0f - c2.y, 1.0f - c2.z, c2.w + opaque);
					shift = 1.0f;
				}

				float x0 = pos.x + container->pos.x + (test(COMPONENT_CENTERED)*(+container->size.x / 2 - size.x / 2));
				float y0 = pos.y + container->pos.y;
				vec4 sPos = vec4();

				// slider button position according to oriantation
				if (vertical)
					sPos = vec4(x0, y0 + float(value) / float(max) * size.y - size.x / 2);
				else
					sPos = vec4(x0 + float(value) / float(max) * size.x - size.y / 2, y0);

				// draw back
				glPushMatrix();
				glColor4f(bColour.x, bColour.y, bColour.z, bColour.w + opaque);
				if (vertical)
				{
					glTranslatef(x0 + size.x / 2 - 1, y0, 0.0f);
					glScalef(3, size.y, 0.0f);
				}
				else
				{
					glTranslatef(x0, y0 + size.y / 2 - 1, 0.0f);
					glScalef(size.x, 3, 0.0f);
				}
				fillRect();
				glPopMatrix();

				// draw border for button
				if (test(COMPONENT_BORDER) && isDown)
				{
					glPushMatrix();
					glColor4f(c2.x, c2.y, c2.z, c2.w + !isDown);
					glTranslatef(sPos.x - 1, sPos.y - 1, 0.0f);
					if (vertical)
						glScalef(size.x + 2, size.x + 2, 0.0f);
					else
						glScalef(size.y + 2, size.y + 2, 0.0f);
					if (test(COMPONENT_ROUNDED))
						fillRRect(0.9f, 10);
					else
						fillRect();
					glPopMatrix();
				}

				// draw button inside
				glPushMatrix();
				glColor4f(c1.x, c1.y, c1.z, c1.w + isDown);
				glTranslatef(sPos.x, sPos.y, 0.0f);
				if (vertical)
					glScalef(size.x, size.x, 0.0f);
				else
					glScalef(size.y, size.y, 0.0f);
				if (test(COMPONENT_ROUNDED))
					fillRRect(0.9f, 10);
				else
					fillRect();
				glPopMatrix();

			}
		}
	}
	void Slider::handleEvents()
	{
		if (container != nullptr)
		{
			if (test(COMPONENT_VISIBLE) && test(COMPONENT_ENABLED) && container->isInside() && !gui->moving)
			{
				bool oldB = test(MOUSE_DOWN);

				if (GetAsyncKeyState(VK_LBUTTON))
					flags |= MOUSE_DOWN;
				else
					flags &= ~MOUSE_DOWN;

				if (isInsideEx() && (gui->id_held == id || gui->id_held == -1) || alreadyDown)
				{
					if (test(MOUSE_DOWN))
					{
						alreadyDown = true;
						gui->id_held = id;

						int oldV = value;
						if (vertical)
							setValue((gui->mousePos.y - pos.y - container->pos.y) / size.y * max);
						else
							setValue((gui->mousePos.x - pos.x - container->pos.x) / size.x * max);

						if (!!(direction = value - oldV) && moveEvent != nullptr)
							moveEvent();
					}
					else
						alreadyDown = false;
				}

				if (gui->id_held == id)
					flags |= COMPONENT_PRESSED;
				else
					flags &= ~COMPONENT_PRESSED;
			}
		}
	}
	bool Slider::isInsideEx()
	{
		// inside button
		bool inside = 0;
		if (container != nullptr)
		{
			// pos of slider button according to oriantation
			vec4 sPos = vec4();
			POINT mp;
			GetCursorPos(&mp);
			mp.x -= glutGet(GLUT_WINDOW_X);
			mp.y -= glutGet(GLUT_WINDOW_Y);

			if (vertical)
			{
				sPos = vec4(pos.x + container->pos.x, pos.y + container->pos.y + value / max * size.y);
				inside = mp.x >= sPos.x && mp.x <= sPos.x + size.x && mp.y >= sPos.y && mp.y <= sPos.y + size.x;
			}
			else
			{
				sPos = vec4(pos.x + container->pos.x + value / max * size.x, pos.y + container->pos.y);
				inside = mp.x >= sPos.x && mp.x < sPos.x + size.y && mp.y >= sPos.y && mp.y <= sPos.y + size.y;
			}
		}
		// inside holder or button
		return isInside() || inside;
	}

	// label
	void Label::handleEvents()
	{
		if (test(COMPONENT_ENABLED | COMPONENT_VISIBLE))
		{

		}
	}
	void Label::setTextMode(int n)
	{
		textMode = n;
	}
	void Label::draw()
	{
		if (container != nullptr)
		{
			if (test(COMPONENT_VISIBLE))
			{
				bool opaque = test(COMPONENT_OPAQUE);

				float x0 = pos.x + container->pos.x + (test(COMPONENT_CENTERED)*(+container->size.x / 2 - size.x / 2));
				float y0 = pos.y + container->pos.y;

				// draw border for back
				if (test(COMPONENT_BORDER))
				{
					glPushMatrix();
					glColor4f(fColour.x, fColour.y, fColour.z, fColour.w + opaque);
					glTranslatef(x0 - 1, y0 - 1, 0.0f);
					glScalef(size.x + 2, size.y + 2, 0.0f);
					if (test(COMPONENT_ROUNDED))
						fillRRect(0.9f, 10);
					else
						fillRect();
					glPopMatrix();
				}

				// draw back
				glPushMatrix();
				glColor4f(bColour.x, bColour.y, bColour.z, bColour.w + opaque);
				glTranslatef(x0, y0, 0.0f);
				glScalef(size.x, size.y, 0.0f);
				if (test(COMPONENT_ROUNDED))
					fillRRect(0.9f, 10);
				else
				{
					if (test(COMPONENT_VERTICAL_GRADIENT))
						vGRect(fColour, bColour);
					else
						fillRect();
				}
				glPopMatrix();

				// draw text
				float l = 0.0f;
				float l2 = 0.0f;
				int tm = 0;
				for (int i = startCursor; i < text.size(); i++, tm++)
				{
					// if wider than label
					if (l >= size.x - glutBitmapWidth(FONT, text[i]) || text[i] == '\n')
					{
						// if new line, remove newline charcter
						if (text[i] == '\n')
						{
							tm = 0;
							i++;
						}
						l = 0.0f;
						l2 += 10;
						// if higher than label
						if (l2 >= size.y - 10)
							break;
					}
					// prints green for first textmode characters
					if (tm < textMode)
						glColor3f(0.2, 0.8, 0.2);
					else
						glColor3f(textColour.x, textColour.y, textColour.z);

					drawTextX(x0 + (test(COMPONENT_CENTERED)*(+container->size.x / 2 - size.x / 2)) + l, y0 + l2 + 10.0f, text[i], FONT);
					l += glutBitmapWidth(FONT, text[i]) + 0.1f;
				}
			}
		}

	}
	void Label::setText(string _text)
	{
		// find the new start cursor
		float l = 0.0f;
		float l2 = 0.0f;
		int i = _text.size() - 1;
		for (; i >= 0; i--)
		{
			if (l >= size.x - glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, _text[i]) || _text[i] == '\n')
			{
				if (l2 >= size.y - 10)
					break;
				l = 0.0f;
				l2 += 10;
				if (_text[i] == '\n')
					i--;
			}
			l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, _text[i]) + 0.1f;
		}
		startCursor = i + 1;
		text = _text;

		// find amount of lines
		l = 0.0f;
		i = 0;
		int _lines = 0;
		for (; i < text.size(); ++i)
		{
			if (l >= size.x - glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, _text[i]) || _text[i] == '\n')
			{
				_lines++;
				l = 0.0f;
				if (_text[i] == '\n')
					i++;
			}
			l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, _text[i]) + 0.1f;
		}
		lines = _lines;
		//cout << "Lines: " << lines << endl;
	}
	void Label::appendText(string _text)
	{
		setText(text + _text);
	}
	void Label::scrollUp()
	{
		if (startCursor != 0)
		{
			// find amount of charactes that fit on the line above the view and subtract amount from startCursor
			float l = 0.0f;
			float l2 = 0.0f;
			int i = startCursor - 1;
			i -= text[i] == '\n';
			for (; i >= 0; i--)
			{
				if (l >= size.x - glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[i]) || text[i] == '\n')
					break;
				l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[i]) + 0.1f;
			}
			startCursor = i + 1;
		}
	}
	void Label::scrollDown()
	{
		// find how many characters on the line below
		float l = 0.0f;
		float l2 = 0.0f;
		int i = startCursor;
		for (; i < text.size(); i++)
		{
			if (l >= size.x - glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[i]) || text[i] == '\n')
				break;
			l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[i]) + 0.1f;
		}

		// check that the bottom line of the text is the bottom line
		l = 0.0f;
		l2 = 0.0f;
		int j = startCursor;
		for (; j < text.size(); j++)
		{
			if (l >= size.x - glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[j]) || text[j] == '\n')
			{
				i += text[i] == '\n';
				l = 0.0f;
				l2 += 10;
				if (l2 >= size.y - 10)
					break;
			}
			l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[j]) + 0.1f;
		}

		if (l2 >= size.y - 10)
			startCursor = i;
	}

	// text field
	int	TextField::findCursorPos()
	{
		POINT mp;
		GetCursorPos(&mp);
		mp.x -= glutGet(GLUT_WINDOW_X);
		mp.y -= glutGet(GLUT_WINDOW_Y);
		// finds what cursor position the mouse is at on the textfield
		if (mp.x >= pos.x + container->pos.x
			&& mp.x <= pos.x + size.x + container->pos.x
			&& mp.y >= pos.y + container->pos.y
			&& mp.y <= pos.y + size.y + container->pos.y)
		{
			float l = 0;
			for (int i = startPos; i < text.size(); i++)
			{
				l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[i]);
				if (pos.x + container->pos.x + l >= mp.x)
					return i;
			}
		}
		return text.length();
	}
	void TextField::handleEvents()
	{
		if (container != nullptr)
		{
			if (test(COMPONENT_VISIBLE) && test(COMPONENT_ENABLED))
			{
				bool oldB = test(MOUSE_DOWN);

				if (GetAsyncKeyState(VK_LBUTTON))
					flags |= MOUSE_DOWN;
				else
					flags &= ~MOUSE_DOWN;

				if (container->isInside())
				{
					if (isInside())
						glutSetCursor(GLUT_CURSOR_TEXT);
					else
						glutSetCursor(GLUT_CURSOR_INHERIT);

					if (isInside() && (gui->id_held == id || gui->id_held == -1))
					{
						// if clicked
						if (!test(MOUSE_DOWN) && oldB)
						{
							if (!test(COMPONENT_FOCUSED))
								flags |= COMPONENT_FOCUSED;
							if (test(MOUSE_DOWN))
								gui->id_held = id;
							cursorPos = findCursorPos();
						}
					}
					else if (test(MOUSE_DOWN))
						flags &= ~COMPONENT_FOCUSED;
				}

				// clamp positions
				cursorPos = cursorPos < 0 ? 0 : cursorPos;
				cursorPos = cursorPos > text.size() ? text.size() : cursorPos;
				startPos = startPos > cursorPos ? cursorPos : startPos;

				if (test(COMPONENT_FOCUSED))
				{
					// test the keydown event
					/// ADD DELETE HELD DOWN EVENT
					if (!(keyboard->keyDownEvent(VK_LBUTTON) || keyboard->keyDownEvent(VK_RBUTTON)))
					{
						if (keyboard->keyDownEvent(VK_BACK))
						{
							if (text.size() > 0 && cursorPos != 0)
							{
								text.erase(cursorPos - 1, 1);
								cursorPos--;
							}
						}
						else if (keyboard->keyDownEvent(VK_LEFT))
						{
							cursorPos--;
						}
						else if (keyboard->keyDownEvent(VK_RIGHT))
						{
							cursorPos++;
						}
						else if (keyboard->keyDownEvent(VK_DELETE))
						{
							text.erase(cursorPos, 1);
						}
						else if (keyboard->keyDownEvent(VK_RETURN) && downEvent != nullptr)
						{
							downEvent();
						}
						else if (keyboard->keyDownEvent(VK_TAB) && keyboard->keyDownEvent(VK_SHIFT))
						{
							gui->updateFocus(this, tabBack);
						}
						else if (keyboard->keyDownEvent(VK_TAB))
						{
							gui->updateFocus(this, tabTo);
						}

						else if (keyboard->keyDownEvent(VK_UP) || keyboard->keyDownEvent(VK_DOWN))
						{
						}
						else /*if an alpha numeric*/
						{
							if (keyboard->getChar() != -1)
							{
								string t;
								t += keyboard->getChar();
								text.insert(cursorPos, t);
								cursorPos++;
							}

						}
						// clamp cursor pos
						cursorPos = cursorPos < 0 ? 0 : cursorPos;
						cursorPos = cursorPos > text.size() ? text.size() : cursorPos;

						if (startPos == cursorPos)
						{
							startPos -= 5;
							startPos = startPos < 0 ? 0 : startPos;
						}

						// calculate new startpos
						if (!test(COMPONENT_MULTILINE))
						{
							float l = 0.0f;
							for (int i = startPos; i < cursorPos; i++)
							{
								l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[i]);
								if (l > size.x)
									startPos++;
							}
						}
						else
						{
							float l = 0.0f;
							float l2 = 0.0f;
							for (int i = startPos; i < cursorPos; i++)
							{
								if (l >= size.x)
								{
									l = 0.0f;
									if (l2 >= size.y - 10)
										startPos++;
									l2 += 10;
								}
								l += glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, text[i]);
							}
						}

						// clamp
						cursorPos = cursorPos < 0 ? 0 : cursorPos;
						cursorPos = cursorPos > text.size() ? text.size() : cursorPos;
						startPos = startPos > cursorPos ? cursorPos : startPos;
					}

				}

			}

		}
	}
	void TextField::draw()
	{

		if (container != nullptr)
		{
			if (test(COMPONENT_VISIBLE))
			{
				bool opaque = test(COMPONENT_OPAQUE);

				float x0 = pos.x + container->pos.x + (test(COMPONENT_CENTERED)*(+container->size.x / 2 - size.x / 2));
				float y0 = pos.y + container->pos.y;

				// border
				if (test(COMPONENT_BORDER))
				{
					glPushMatrix();
					glColor4f(fColour.x, fColour.y, fColour.z, fColour.w + opaque);
					glColor4f(fColour.x, fColour.y, fColour.z, fColour.w + opaque);
					glTranslatef(x0 - 1, y0 - 1, 0.0f);
					glScalef(size.x + 2, size.y + 2, 0.0f);
					if (test(COMPONENT_ROUNDED))
						fillRRect(0.9f, 10);
					else
						fillRect();
					glPopMatrix();
				}

				// bg
				glPushMatrix();
				glColor4f(bColour.x, bColour.y, bColour.z, bColour.w + opaque);
				glTranslatef(x0, y0, 0.0f);
				glScalef(size.x, size.y, 0.0f);
				if (test(COMPONENT_ROUNDED))
					fillRRect(0.9f, 10);
				else
				{
					if (test(COMPONENT_VERTICAL_GRADIENT))
						vGRect(fColour, bColour);
					else
						fillRect();
				}
				glPopMatrix();

				// text
				if (!test(COMPONENT_MULTILINE))
				{
					float l = 0.0f;
					endPos = 0;
					for (int i = startPos; i < text.size(); i++, endPos++)
					{
						if (l >= size.x - glutBitmapWidth(FONT, test(COMPONENT_CENSORED) ? '*' : text[i]))
							break;
						l += glutBitmapWidth(FONT, test(COMPONENT_CENSORED) ? '*' : text[i]);
					}
					glColor3f(textColour.x, textColour.y, textColour.z);
					float x1 = x0;
					float y1 = y0 + (size.y - 10) / 2.0f;
					drawTextX(x1, y1 + 10.0f, test(COMPONENT_CENSORED) ? repeatChar('*', endPos) : text.substr(startPos, endPos), FONT);
					//glPopMatrix();

					if (test(COMPONENT_FOCUSED) && int(gui->blink) % 2)
					{
						//glPushMatrix();
						glColor3f(textColour.x, textColour.y, textColour.z);
						l = 0.0f;
						for (int i = startPos; i < cursorPos; i++)
						{
							if (l >= size.x - glutBitmapWidth(FONT, test(COMPONENT_CENSORED) ? '*' : text[i]))
								break;
							l += glutBitmapWidth(FONT, test(COMPONENT_CENSORED) ? '*' : text[i]);
						}
						x1 = x1 = x0 - 1.0f + l;
						y1 = y0 + (size.y - 10 - 3.0f) / 2.0f;
						drawTextX(x1, y1 - 7.0f + 10.0f, "|", FONT);
						//glPopMatrix();
					}
				}
				else
				{
					float x1 = x0;
					float y1 = y0 + 3.0f;
					float l = 0.0f;
					float l2 = 0.0f;
					for (int i = startPos; i < text.size(); i++)
					{
						if (l >= size.x - glutBitmapWidth(FONT, text[i]))
						{
							l = 0.0f;
							if (l2 >= size.y - 10)
								break;
							l2 += 10;
						}
						glColor3f(textColour.x, textColour.y, textColour.z);
						drawTextX(x1 + l, y1 + l2 + 10.0f, text.substr(i, 1), FONT);
						//glPopMatrix();
						l += glutBitmapWidth(FONT, text[i]) + 0.1f;
					}

					if (test(COMPONENT_FOCUSED) && int(gui->blink) % 2)
					{
						//glPushMatrix();
						glColor3f(textColour.x, textColour.y, textColour.z);
						x1 = x0 - 1.0f;
						y1 = y0 - 5.0f;
						l = 0.0f;
						l2 = 0.0f;
						for (int i = startPos; i < cursorPos; i++)
						{
							if (l >= size.x - glutBitmapWidth(FONT, text[i]) || text[i] == '\n')
							{
								l = 0.0f;
								if (l2 >= size.y - 10)
									break;
								l2 += 10;
							}
							l += glutBitmapWidth(FONT, text[i]);
						}
						x1 += l;
						y1 += l2;
						drawTextX(x1, y1 + 10.0f, "|", FONT);
						//glPopMatrix();
					}
				}
			}
		}
	}
	
	// button
	void Button::handleEvents()
	{
		if (container != nullptr)
		{
			if (test(COMPONENT_VISIBLE) && test(COMPONENT_ENABLED) && container->isInside() && !gui->moving)
			{
				bool oldB = test(MOUSE_DOWN);

				if (GetAsyncKeyState(VK_LBUTTON))
					flags |= MOUSE_DOWN;
				else
					flags &= ~MOUSE_DOWN;

				if (isInside() && (gui->id_held == id || gui->id_held == -1))
				{
					// if clicked
					if (!test(MOUSE_DOWN) && oldB)
					{
						if (!test(COMPONENT_TOGGLED))
						{
							if (downEvent != nullptr) downEvent();
						}
						else
							if (upEvent != nullptr) upEvent();
						if (test(COMPONENT_TOGGLABLE)) flags ^= COMPONENT_TOGGLED;
					}
					if (test(MOUSE_DOWN))
						gui->id_held = id;
				}

				if (gui->id_held == id)
					flags |= COMPONENT_PRESSED;
				else
					flags &= ~COMPONENT_PRESSED;
			}
		}
	}
	void Button::draw()
	{
		if (container != nullptr)
		{
			if (test(COMPONENT_VISIBLE))
			{
				bool opaque = test(COMPONENT_OPAQUE);
				bool isDown = test(COMPONENT_PRESSED) || test(COMPONENT_TOGGLED);
				vec4 c1 = bColour;
				vec4 c2 = fColour;
				vec4 c3 = textColour;
				float shift = 0.0f;
				if (isDown)
				{
					c1 = vec4(1.0f - c1.x, 1.0f - c1.y, 1.0f - c1.z, c1.w + opaque);
					c2 = vec4(1.0f - c2.x, 1.0f - c2.y, 1.0f - c2.z, c2.w + opaque);
					c3 = vec4(1.0f - c3.x, 1.0f - c3.y, 1.0f - c3.z, c3.w + opaque);
					shift = 1.0f;
				}

				// border
				float x0 = pos.x + container->pos.x + (test(COMPONENT_CENTERED)*(+container->size.x / 2 - size.x / 2));
				float y0 = pos.y + container->pos.y;

				if (test(COMPONENT_BORDER) && isDown)
				{
					glPushMatrix();
					glColor4f(c2.x, c2.y, c2.z, c2.w + !isDown);
					glTranslatef(x0 - 1, y0 - 1, 0.0f);
					glScalef(size.x + 2, size.y + 2, 0.0f);
					if (test(COMPONENT_ROUNDED))
						fillRRect(0.9f, 10);
					else
						fillRect();
					glPopMatrix();
				}

				// middle
				glPushMatrix();
				glColor4f(c1.x, c1.y, c1.z, c1.w + isDown);
				glTranslatef(x0, y0, 0.0f);
				glScalef(size.x, size.y, 0.0f);
				if (test(COMPONENT_ROUNDED))
					fillRRect(0.9f, 10);
				else if (test(COMPONENT_VERTICAL_GRADIENT))
					vGRect(c1, c2);
				else if (test(COMPONENT_HORIZONTAL_GRADIENT))
					hGRect(c1, c2);
				else
					fillRect();
				glPopMatrix();

				// text
				string print_str = isDown ? downText : text;
				float l = 0;
				for (int i = 0; i < print_str.length(); i++)
					l += glutBitmapWidth(FONT, print_str[i]);
				glPushMatrix();
				glColor3f(c3.x, c3.y, c3.z);
				float x1 = x0 + (size.x - l - 2.0f) / 2.0f + shift;
				float y1 = y0 + (size.y - 10) / 2.0f + shift;
				drawTextX(x1, y1 + 10.0f, print_str, FONT);
				glPopMatrix();
			}
		}
	}
	void Button::setText(string txt)
	{
		text = txt;
		downText = txt;
	}
	void Button::setText(string upTxt, string downTxt)
	{
		text = upTxt;
		downText = downTxt;
	}

	// gui manager
	void GUIManager::handleEvents()
	{
		keyboard->run();
		blink += 0.02f;
		amountOn = 0;
		for (int i = 0; i < components.size(); ++i)
		{
			components[i]->handleEvents();
			amountOn += components[i]->isInside() && components[i]->test(COMPONENT_VISIBLE);
		}
		mouseOnGUI = amountOn > 0;
		id_held = GetAsyncKeyState(VK_LBUTTON) ? id_held : -1;
		updateMousePosition();
		windowMove();
		focusChange();
		/*cDown = 0;
		for (int i = 0; i < components.size(); ++i)
		cDown += components[i]->isInside();
		mouseOnGUI = cDown > 0;*/
	}
	void GUIManager::draw()
	{
		for (int i = components.size() - 1; i >= 0; --i)
			components[i]->draw();
	}
	void GUIManager::add(Component *c)
	{
		c->container = background;
		c->id = ids++;
		c->gui = this;
		components.push_back(c);
	}
	void GUIManager::updateMousePosition()
	{
		GetCursorPos(&cursorPos);
		mousePos = vec4(cursorPos.x - glutGet(GLUT_WINDOW_X), cursorPos.y - glutGet(GLUT_WINDOW_Y));
	}
	void GUIManager::windowMove()
	{
		vec4 windowSize = vec4(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
		if (GetAsyncKeyState(VK_LBUTTON))
		{
			// if inside the panel
			if (mousePos.x > 0 && mousePos.x < windowSize.x && mousePos.y > 0 && mousePos.y < 0)
			{
				moving = true;
				if (mouseOffset.x == 0 && mouseOffset.y == 0)
					mouseOffset = mousePos;
			}
		}
		else /*if released*/
		{
			moving = false;
			mouseOffset = vec4();
		}

		if (moving && id_held == -1)
			glutPositionWindow(cursorPos.x - mouseOffset.x, cursorPos.y - mouseOffset.y);
	}
	void GUIManager::focusChange()
	{
		if (from != nullptr && to != nullptr)
		{
			from->flags &= ~COMPONENT_FOCUSED;
			to->flags |= COMPONENT_FOCUSED;
		}
		from = nullptr;
		to = nullptr;
	}
	void GUIManager::updateFocus(Component * _from, Component * _to)
	{
		from = _from;
		to = _to;
	}

	// OTHER GUI FUNCTIONS //
	void removeWindowBorder()
	{
		HWND hwnd = FindWindowA(NULL, gui->windowName.c_str());
		DWORD style = GetWindowLong(hwnd, GWL_STYLE);
		style &= ~WS_OVERLAPPEDWINDOW;
		SetWindowLong(hwnd, GWL_STYLE, style);
	}
	void toggleGlutWindowMaximizeBox()
	{
		HWND hwnd = FindWindowA(NULL, gui->windowName.c_str());
		DWORD style = GetWindowLong(hwnd, GWL_STYLE);
		style = WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX;
		SetWindowLong(hwnd, GWL_STYLE, style);
	}
	void minimise()
	{
		ShowWindow(FindWindowA(NULL, gui->windowName.c_str()), SW_MINIMIZE);
	}
	void killApp()
	{
		exit(EXIT_SUCCESS);
	}

	void handlerThread()
	{
		while (1)
		{
			if (GetForegroundWindow() == win)
				gui->handleEvents();
			// reduce load
			this_thread::sleep_for(std::chrono::milliseconds(10));
		}
	}
	void detachEventHandler()
	{
		if (!gui->locked)
		{
			gui->locked = true;
			thread t1(handlerThread);
			t1.detach();
		}
		else
		{
			cerr << "ERROR: Handler thread already running." << endl;
			throw exception("ERROR: Handler thread already running.");
		}
	}

}

#endif