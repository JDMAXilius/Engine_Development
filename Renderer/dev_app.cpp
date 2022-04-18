#include "dev_app.h"
#include "math_types.h"
#include "debug_renderer.h"
#include <iostream>
#include <ctime>
#include <random>
#include <math.h>
#include <bitset>
#include <Windows.h>
#include <DirectXMath.h>
#include "frustum_culling.h"

#define PI 3.14159265

//TODO include debug_renderer.h and pools.h and anything else you might need here
#include "pools.h"
using namespace std;
using namespace DirectX;

namespace end
{
	std::array<aabb_t, 12> boxesCollision;
	frustum_t frustum;
	points frustumCorners;

	float3 flot3;
	
	struct particle
	{
		float3 pos;
		float3 prev_pos;
		float4 color;
		//could add lifetime if you’d like
		float3 Velocity;
		double lifetime;
	};
	particle pt;
	particle pt2;
	particle pt3;
	particle pt4;

	//srand();
	//float x = (float)rand() / ((float)RAND_MAX / 1);
	pool_t<particle, 2024*4> shared_pool;
	sorted_pool_t<particle, 1024> sorted_pool;
	struct emitter
	{
		float3 spawn_pos;
		//float3 spawn_prevpos;
		float4 spawn_color;
		//float3 Velocity;
		// indices into the shared_pool 
		sorted_pool_t<int16_t, 1024> indices;
	};
	emitter fp;
	emitter fp2;
	emitter fp3;
	emitter fp4;

	//Particle behaviour
	float3 Gravity = float3(0, -9.8f, 0);
	float3 Pos;
	
	//VK_DOWN
	//VK_RIGHT

	double delta_time = 0.0, times = 0.0;
	int x=0, y = 0, z = 0;

	//Matrix Behaviors
	float4x4  Look_At, Turn_To, Player;
	float4 X_Axis, Y_Axis, Z_Axis;
	float4 R = { 1.0f, 0,0, 1 };
	float4 G = { 0, 1, 0, 1 };
	float4 B = { 0, 0, 1, 1 };
	//XMMATRIX LookAt, TurnTo;

	float4x4 identityD()
	{
		float4x4 t = { float4(0,0,0,0) };
		t[0].x = 1;
		t[1].y = 1;
		t[2].z = 1;
		t[3].w = 1;
		return t;
	}
	float4x4 identity =
	{
		float4(1, 0, 0, 0),
		float4(0, 1, 0, 0),
		float4(0, 0, 1, 0),
		float4(0, 0, 0, 1)
	};
	 
	double dev_app_t::get_delta_time()const
	{
		return delta_time;
	}

	
	float RandomFloat(float a, float b) {
		float random = ((float)rand()) / (float)RAND_MAX;
		float diff = b - a;
		float r = random * diff;
		return a + r;
	}
	dev_app_t::dev_app_t()
	{
		Player = identity;
		Look_At = identity;
		Turn_To = identity;

		//Player[3].xyz = float3(0, 0, 0);
		Look_At[3].xyz = float3(-5, 2, -5);
		Turn_To[3].xyz = float3(5, 2, 3);

		//std::bitset;
		//GetAsyncKeyState(…),
		int i = 0;
		while (i < boxesCollision.size())
		{
			boxesCollision[i].center = float3(RandomFloat(-14, 14), RandomFloat(1, 2), RandomFloat(-14, 14));
			boxesCollision[i].extents = float3(RandomFloat(0.2f, 2.2), RandomFloat(0.2f, 2.2), RandomFloat(0.2f, 2.2));
			i++;
		}

		std::cout << "Log whatever you need here.\n"; // Don’t forget to include <iostream>
	}
	
	double calc_delta_time() 
	{
		static std::chrono::time_point<std::chrono::high_resolution_clock> last_time = std::chrono::high_resolution_clock::now();

		std::chrono::time_point<std::chrono::high_resolution_clock> new_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_seconds = new_time - last_time;
		last_time = new_time;

		return min(1.0 / 15.0, elapsed_seconds.count());
	}

	void grid(double delta)
	{
		z = 0;
		x = 0;
		for (int i = 0; i < 11; i++)
		{
			//end::debug_renderer::add_line(float3(-10, 0, z), float3(10, 0, z), float4(0.1f, 1, 0.1f, 1));
			end::debug_renderer::add_line(float3(-10, 0, z), float3(10, 0, z), float4((sin(delta) / (PI / 180)), (cos(delta) / (PI / 180)), sin(delta) / 4 + 0.5, 1));
			z++;
			//end::debug_renderer::add_line(float3(x, 0, -10), float3(x, 0, 10), float4(0.1f, 1, 0.1f, 1));
			end::debug_renderer::add_line(float3(x, 0, -10), float3(x, 0, 10), float4((sin(delta) / (PI / 180)), (cos(delta ) / (PI / 180)), sin(delta) / 4 + 0.5, 1));
			x++;
		}
		z = 0;
		x = 0;
		for (int i = 0; i < 11; i++)
		{
			//end::debug_renderer::add_line(float3(-10, 0, z), float3(10, 0, z), float4(0.1f, 1, 0.1f, 1));
			end::debug_renderer::add_line(float3(-10, 0, z), float3(10, 0, z), float4((sin(delta) / (PI / 180)), (cos(delta) / (PI / 180)), sin(delta) / 4 + 0.5, 1));
			z--;
			//end::debug_renderer::add_line(float3(x, 0, -10), float3(x, 0, 10), float4(0.1f, 1, 0.1f, 1));
			end::debug_renderer::add_line(float3(x, 0, -10), float3(x, 0, 10), float4((sin(delta) / (PI / 180)), (cos(delta) / (PI / 180)), sin(delta) / 4 + 0.5, 1));
			x--;
		}
	}

	void sortedEmiter(particle pt, int number)
	{
		for (size_t i = 0; i < number; i++)
		{
			int index = sorted_pool.alloc();
			if (index != -1)
			{
				// first initialize
				sorted_pool[index].pos = pt.pos;
				//sorted_pool[index].prev_pos = pt.prev_pos;
				sorted_pool[index].color = pt.color;
				sorted_pool[index].Velocity = float3(RandomFloat(-1, 1), RandomFloat(8, 12), RandomFloat(-1, 1));
				//sorted_pool[index] = pt;
			}
		}

		for (size_t i = 0; i < sorted_pool.size(); i++)
		{
			sorted_pool[i].lifetime -= delta_time;
			// then update
			if (sorted_pool[i].pos.y >= 0)
			{
				sorted_pool[i].prev_pos = sorted_pool[i].pos;
				sorted_pool[i].pos += (sorted_pool[i].Velocity * delta_time) + (Gravity * 0.5f * pow(delta_time, 2));
				sorted_pool[i].Velocity += Gravity * delta_time;
			}
			else if(sorted_pool[i].lifetime <= 0)
			{
				sorted_pool.free(i);
				i--;
			}
		}

		for (size_t i = 0; i < sorted_pool.size(); i++)
		{
			// lastly render
			end::debug_renderer::add_line(sorted_pool[i].prev_pos, sorted_pool[i].pos, sorted_pool[i].color);
		}

	}

	void freeEmiter(emitter &pt, int number)
	{
		for (size_t i = 0; i < number; i++)
		{
			int16_t fpIndex = shared_pool.alloc();
			if (fpIndex != -1)
			{
				
				int16_t spIndex = pt.indices.alloc();
				if (spIndex != -1)
				{
					// first initialize
					shared_pool[fpIndex].pos = pt.spawn_pos;
					//shared_pool[index].prev_pos = pt.spawn_prevpos;
					shared_pool[fpIndex].color = pt.spawn_color;
					shared_pool[fpIndex].Velocity = float3(RandomFloat(-1, 1), RandomFloat(8, 12), RandomFloat(-1, 1));
					//sorted_pool[index] = pt;
					pt.indices[spIndex] = fpIndex;
				}
				else
				{
					shared_pool.free(fpIndex);
				}
			}
		}

		for (size_t i = 0; i < pt.indices.size(); i++)
		{
			shared_pool[i].lifetime -= delta_time;
			// then update
			if (shared_pool[pt.indices[i]].pos.y >= 0)
			{
				shared_pool[pt.indices[i]].prev_pos = shared_pool[pt.indices[i]].pos;
				shared_pool[pt.indices[i]].pos += (shared_pool[pt.indices[i]].Velocity * delta_time) + (Gravity * 0.5f * pow(delta_time, 2));
				shared_pool[pt.indices[i]].Velocity += Gravity * delta_time;
			}
			else if (sorted_pool[i].lifetime <= 0)
			{
				shared_pool.free(pt.indices[i]);
				pt.indices.free(i);
				i--;
			}
		}

		for (size_t i = 0; i < pt.indices.size(); i++)
		{
			// lastly render
			end::debug_renderer::add_line(shared_pool[pt.indices[i]].prev_pos, shared_pool[pt.indices[i]].pos, shared_pool[pt.indices[i]].color);
		}

	}
	//
	float3 position_init;
	float3 position_end;
	float4 color;
	float3 dir;
	float scale;
	
	//bitset<128>bits;
	bitset<128>mbits;

	POINT mPosition;
	POINT prevmPosition;

	void ViewLines(float4x4 p, float3 axis, float4 color, int num )
	{

		//Vector translation
		position_init = p[3].xyz;
		position_end;
		//color = { axis.x, axis.y, axis.z, 1 };
		dir = {p[0].x,p[1].y,p[2].z };
		scale = 1;
		//position_init += (dir * scale);
		position_end = position_init + (axis * num);
		end::debug_renderer::add_line(position_init, position_end, color);
	}
	float4x4 matrixMultMatrix4x4(float4x4 m1, float4x4 m2)
	{
		float4x4 mat = identity;;
		
		int i = 0;
		while (i <= 3)
		{
			mat[0][i] = (m1[0][0] * m2[0][i]) + (m1[0][1] * m2[1][i]) + (m1[0][2] * m2[2][i]) + (m1[0][3] * m2[3][i]);

			mat[1][i] = (m1[1][0] * m2[0][i]) + (m1[1][1] * m2[1][i]) + (m1[1][2] * m2[2][i]) + (m1[1][3] * m2[3][i]);

			mat[2][i] = (m1[2][0] * m2[0][i]) + (m1[2][1] * m2[1][i]) + (m1[2][2] * m2[2][i]) + (m1[2][3] * m2[3][i]);

			mat[3][i] = (m1[3][0] * m2[0][i]) + (m1[3][1] * m2[1][i]) + (m1[3][2] * m2[2][i]) + (m1[3][3] * m2[3][i]);
			i++;
		};

		return mat;
	}
	float4x4 rotateY(float4x4 p , double d)
	{
		d = d * 3.14f / 180.0f;
		float4x4 t = identity;
		t[0].x = cos(d);
		t[2].x = sin(d);
		t[0].z = -sin(d);
		t[2].z = cos(d);
		t[1].y = 1;
		t[3].w = 1;
		
		return matrixMultMatrix4x4(t, p);
	}
	float dot(float3 a, float3 b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	float3 cross(float3 a, float3 b)
	{
		float3 vec = { (a.y * b.z - a.z * b.y), (a.z * b.x - a.x * b.z), (a.x * b.y - a.y * b.x)};
		return vec;
	}
	float vlength(float3 vec)
	{
		return sqrtf(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
	}
	float3 NormalizeVector(float3 vec) {
		float magnitude = vlength(vec);
		if (magnitude != 0)
		{
			vec.x /= magnitude;
			vec.y /= magnitude;
			vec.z /= magnitude;
		}
		else
			vec.x = vec.y = vec.z = 0.0f;

		return vec;
	}

	//float4x4 turnTo(float4x4 L, float3 T)
	//{
	//	//float4x4 v = (L[3].w - T[3].w).normalize();
	//	float3 point = (T - L[3].xyz);
	//	//NormalizeVector(point);
	//	point = flot3.normalize(point);
	//	//XVECTOR 
	//	float dotY = flot3.dot(point, Turn_To[0].xyz);
	//	
	//	//float4x4 RY = rotateY(dotY * delta_time);
	//	//float4x4 m = { float4(x.x,x.y,x.z,0), float4(y.x,y.y,y.z,0), float4(z.x,z.y,z.z,0), L[3] };
	//	//L = RY(L);
	//	//dotX = dot(v, L[1].xyz);
	//	//RX = matrixX(dotX + delta_time);
	//	//L = RX(L);
	//	//return orthonormalize(L);
	//}

	float4x4 orthonormalize(float4x4 L)
	{
		//XMVECTOR wUP = XMVectorSet(0, 1, 0, 0);
		//z = L[2];
		//XMVECTOR zD = XMVector3Normalize(Look_At[2].xyz);
		float3 z = flot3.normalize(Look_At[2].xyz);
		float3 wUP = float3(0, 1, 0);
		float3 x = flot3.normalize(flot3.cross(wUP, z));
		float3 y = flot3.normalize(flot3.cross(z, x));
		float4x4 m = {float4(x.x,x.y,x.z,0), float4(y.x,y.y,y.z,0), float4(z.x,z.y,z.z,0), L[3]};
		return m;		//return float4x4{ x,y,z,Look_At[3].xyz };
	}
	XMMATRIX orthonormalizeD(XMMATRIX L)
	{
		XMVECTOR z = XMVector3Normalize(L.r[2]);
		XMVECTOR wUP = XMVectorSet(0, 1, 0, 0);
		XMVECTOR x = XMVector3Normalize(XMVector3Cross(wUP, z));
		XMVECTOR y = XMVector3Normalize(XMVector3Cross(z, x));
		XMMATRIX m = { x, y, z, L.r[3] };
		return m;
	}
	XMMATRIX turnToD(XMMATRIX L, XMVECTOR T)
	{
		//float4x4 v = (L[3].w - T[3].w).normalize();
		XMVECTOR point = XMVector3Normalize(T - L.r[3]);
		
		XMVECTOR dotY = XMVector3Dot(point, L.r[0]);
		XMMATRIX RY = XMMatrixRotationY(XMVectorGetX(dotY) * delta_time);
		
		XMVECTOR dotX = XMVector3Dot(point, L.r[1]);
		XMMATRIX RX = XMMatrixRotationX(XMVectorGetX(dotX) *  -delta_time);
		L = XMMatrixMultiply(RY, L);
		L = XMMatrixMultiply(RX, L);
		return orthonormalizeD(L);
	}
	
	XMMATRIX loockATD(XMMATRIX L, XMVECTOR T)
	{
		////float4x4 v = (L[3].w - T[3].w).normalize();
		//XMVECTOR point = XMVector3Normalize(T - L.r[3]);

		//XMVECTOR dotZ = XMVector3Dot(point, L.r[2]);
		//XMMATRIX RZ = XMMatrixRotationY(XMVectorGetX(dotZ) * delta_time);
		XMVECTOR z = XMVector3Normalize(T - L.r[3]);
		XMVECTOR wUP = XMVectorSet(0, 1, 0, 0);
		XMVECTOR x = XMVector3Normalize(XMVector3Cross(wUP, z));
		XMVECTOR y = XMVector3Normalize(XMVector3Cross(z, x));
		XMMATRIX m = { x, y, z, L.r[3] };

		return m;
	}
	view_t cameraMove(view_t defaultView, int num, int index, bool b)
	{
		XMMATRIX temp = (XMMATRIX&)defaultView.view_mat;
		XMVECTOR direction = temp.r[index];
		if (b)
		{
			temp.r[3] += direction * delta_time * num;
		}
		else
		{
			temp.r[3] -= direction * delta_time * num;
		}
		defaultView.view_mat = (float4x4_a&)temp;
		return defaultView;
	}
	view_t mouseMove(view_t defaultView, int num)
	{
		XMMATRIX temp = (XMMATRIX&)defaultView.view_mat;
		temp.r[3] -= XMVectorSet(0, 0, 0, 1);
		float dy = mPosition.y - prevmPosition.y;
		float dx = mPosition.x - prevmPosition.x;
		temp = XMMatrixMultiply(temp, XMMatrixRotationY(delta_time * num * dx));
		temp = XMMatrixMultiply(XMMatrixRotationX(delta_time * num * dy), temp);
		temp.r[3] = (XMVECTOR&)defaultView.view_mat[3];
		defaultView.view_mat = (float4x4_a&)temp;
		return defaultView;
	}

	void dev_app_t::setMDown(uint8_t number)
	{
		mbits.set(0, true);
	}
	void dev_app_t::setMUp(uint8_t number)
	{
		mbits.set(0, false);
		//bool check = false;
	}
	void dev_app_t::setMPosition(float x, float y)
	{
		mPosition = { (LONG)x, (LONG)y };
	}

	// Calculates the plane of a triangle from three points.
	plane_t calculate_plane(float3 a, float3 b, float3 c)
	{
		float3 second = (float3&)DirectX::XMVector3Normalize((DirectX::XMVECTOR&)(c - a));
		float3 first = (float3&)DirectX::XMVector3Normalize((DirectX::XMVECTOR&)(b - a));
		float3 normal = (float3&)DirectX::XMVector3Normalize(DirectX::XMVector3Cross((DirectX::XMVECTOR&)first, (DirectX::XMVECTOR&)second));
		float offset = (float&)DirectX::XMVector3Dot((DirectX::XMVECTOR&)a, (DirectX::XMVECTOR&)normal);
		plane_t answer = { normal, offset };
		return answer;
	}

	void renderCalc(float3 nearc, float3 farc, points point)
	{
		float3 center = (point[0] + point[2] + point[4] + point[6]) / 4;
		//debug_renderer::add_line(nc, fc, float4(1,1,0,1), float4(1,1,0,1));
		debug_renderer::add_line(farc, farc + frustum[FARe].normal, float4(1, 1, 1, 1), float4(0, 1, 1, 1));
		debug_renderer::add_line(nearc, nearc + frustum[NEARe].normal, float4(1, 1, 1, 1), float4(0, 1, 1, 1));
		debug_renderer::add_line(center, center + frustum[LEFT].normal, float4(1, 1, 1, 1), float4(0, 1, 1, 1));

		center = (point[1] + point[3] + point[5] + point[7]) / 4;
		debug_renderer::add_line(center, center + frustum[RIGHT].normal, float4(1, 1, 1, 1), float4(0, 1, 1, 1));

		center = (point[2] + point[3] + point[6] + point[7]) / 4;
		debug_renderer::add_line(center, center + frustum[BOTTOM].normal, float4(1, 1, 1, 1), float4(0, 1, 1, 1));

		center = (point[0] + point[1] + point[4] + point[5]) / 4;
		debug_renderer::add_line(center, center + frustum[TOP].normal, float4(1, 1, 1, 1), float4(0, 1, 1, 1));
	}

	void pointCalc(float Hfar, float Wfar, float Hnear, float Wnear, float4x4& m, points point, float3 nearc, float3 farc)
	{
		float Wide = 0.75f;
		float Hieght = 0.5f;
		point[FTR] = (farc + m[1].xyz * (Hfar * Hieght)) + (m[0].xyz * (Wfar * Wide));
		point[FBR] = (farc - m[1].xyz * (Hfar * Hieght)) + (m[0].xyz * (Wfar * Wide));
		point[FBL] = (farc - m[1].xyz * (Hfar * Hieght)) - (m[0].xyz * (Wfar * Wide));
		point[FTL] = (farc + m[1].xyz * (Hfar * Hieght)) - (m[0].xyz * (Wfar * Wide));

		point[NTR] = (nearc + m[1].xyz * (Hnear * Hieght)) + (m[0].xyz * (Wnear * Wide));
		point[NBR] = (nearc - m[1].xyz * (Hnear * Hieght)) + (m[0].xyz * (Wnear * Wide));
		point[NBL] = (nearc - m[1].xyz * (Hnear * Hieght)) - (m[0].xyz * (Wnear * Wide));
		point[NTL] = (nearc + m[1].xyz * (Hnear * Hieght)) - (m[0].xyz * (Wnear * Wide));
	}
	void frustumCal(points point, frustum_t& frustum)
	{
		frustum[TOP] = calculate_plane(point[FTR], point[FTL], point[NTL]);
		frustum[RIGHT] = calculate_plane(point[FTR], point[NBR], point[FBR]);
		frustum[LEFT] = calculate_plane(point[NTL], point[FTL], point[NBL]);
		frustum[BOTTOM] = calculate_plane(point[FBL], point[FBR], point[NBR]);
		frustum[FARe] = calculate_plane(point[FTL], point[FTR], point[FBL]);
		frustum[NEARe] = calculate_plane(point[NTR], point[NTL], point[NBL]);
	}

	// Calculates a frustum (6 planes) from the input view parameter.
	void calculate_frustum(frustum_t& frustum, float4x4& m, points& point)
	{
		float fov = DirectX::XMConvertToRadians(25);
		float ar = 1.5;
		float nearDist = 1;
		float farDist = 10;
		
		float3 nearc = m[3].xyz + m[2].xyz * nearDist;
		float3 farc = m[3].xyz + m[2].xyz * farDist;

		float Wide = 0.75f;
		float Hnear = 2 * tan(fov / 2) * nearDist;
		float Hfar = 2 * tan(fov / 2) * farDist;
		float Hieght = 0.5f;
		float Wnear = Hnear * ar;
		float Wfar = Hfar * ar;

		point[FTL] = (farc + m[1].xyz * (Hfar * Hieght)) - (m[0].xyz * (Wfar * Wide));
		//pointCalc(fc, Hfar, Wfar, m, FTL, point);
		point[FTR] = (farc + m[1].xyz * (Hfar * Hieght)) + (m[0].xyz * (Wfar * Wide));
		point[FBL] = (farc - m[1].xyz * (Hfar * Hieght)) - (m[0].xyz * (Wfar * Wide));
		point[FBR] = (farc - m[1].xyz * (Hfar * Hieght)) + (m[0].xyz * (Wfar * Wide));

		point[NTL] = (nearc + m[1].xyz * (Hnear * Hieght)) - (m[0].xyz * (Wnear *Wide));
		point[NTR] = (nearc + m[1].xyz * (Hnear * Hieght)) + (m[0].xyz * (Wnear *Wide));
		point[NBL] = (nearc - m[1].xyz * (Hnear * Hieght)) - (m[0].xyz * (Wnear *Wide));
		point[NBR] = (nearc - m[1].xyz * (Hnear * Hieght)) + (m[0].xyz * (Wnear *Wide));
		
		//pointCalc( Hfar, Wfar, Hnear, Wnear,  m, point, nearc, farc);

		frustumCal(point, frustum);

		renderCalc( nearc, farc, point);
	}

	// Calculates which side of a plane the sphere is on.
	int classify_sphere_to_plane(const sphere_t& sphere, const plane_t& plane)
	{
		float soffset = (float&)DirectX::XMVector3Dot((DirectX::XMVECTOR&)sphere.center, (DirectX::XMVECTOR&)plane.normal);
		
		float balltoplane = soffset - plane.offset;
		
		int STP = (balltoplane > sphere.radius) ? 1 : ((balltoplane < -sphere.radius) ? -1 : 0);
		return STP;
	}
	sphere_t sphere;
	float3 absnormal;
	// Calculates which side of a plane the aabb is on.
	int classify_aabb_to_plane(const aabb_t& aabb, const plane_t& plane)
	{
		absnormal = { abs(plane.normal.x), abs(plane.normal.y), abs(plane.normal.z) };
		float radius = (float&)DirectX::XMVector3Dot((DirectX::XMVECTOR&)aabb.extents, (DirectX::XMVECTOR&)absnormal);
		sphere.radius = radius;
		sphere.center = aabb.center;
		int STP = classify_sphere_to_plane(sphere, plane);
		return STP;
	}

	// Determines if the aabb is inside the frustum.
	bool aabb_to_frustum(const aabb_t& aabb, const frustum_t& frustum)
	{
		int i = 0;
		while (i < 6)	
		{
			int STP = classify_aabb_to_plane(aabb, frustum[i]);
			i++;
			if (STP == -1)
			{
				return false;
			}
		}
		return true;
	}


	view_t dev_app_t::update(view_t defaultView)
	{
		delta_time = calc_delta_time();
		times += delta_time;

		ViewLines(Player, Player[0].xyz, R, 2);
		ViewLines(Player, Player[1].xyz, G, 2);
		ViewLines(Player, Player[2].xyz, B, 2);

		ViewLines(Look_At, Look_At[0].xyz, R, 2);
		ViewLines(Look_At, Look_At[1].xyz, G, 2);
		ViewLines(Look_At, Look_At[2].xyz, B, 2);

		ViewLines(Turn_To, Turn_To[0].xyz, R, 2);
		ViewLines(Turn_To, Turn_To[1].xyz, G, 2);
		ViewLines(Turn_To, Turn_To[2].xyz, B, 2);

		/*if (GetAsyncKeyState('D') & 0x1 && (GetAsyncKeyState('W') & 0x1)) {
			Player = rotateY(Player, 1000 * delta_time);
			Player[3].xyz += Player[2].xyz * 100 * delta_time;
		}
		if ((GetAsyncKeyState('A') & 0x1) && (GetAsyncKeyState('W') & 0x1)) {
			Player = rotateY(Player, -1000 * delta_time);
			Player[3].xyz += Player[2].xyz * 100 * delta_time;
		}
		if (GetAsyncKeyState('D') & 0x1 && (GetAsyncKeyState('S') & 0x1)) {
			Player = rotateY(Player, 1000 * delta_time);
			Player[3].xyz += Player[2].xyz * -100 * delta_time;
		}
		if ((GetAsyncKeyState('A') & 0x1) && (GetAsyncKeyState('S') & 0x1)) {
			Player = rotateY(Player, -1000 * delta_time);
			Player[3].xyz += Player[2].xyz * -100 * delta_time;
		}*/
		if (GetAsyncKeyState(VK_UP) & 0x1) 
		{
			Player[3].xyz += Player[2].xyz * 200 * delta_time;
			//cameraMove( *defaultView);
		}
		if (GetAsyncKeyState(VK_DOWN) & 0x1) 
		{
			Player[3].xyz += Player[2].xyz * -200 * delta_time;
			
		}
		if (GetAsyncKeyState(VK_RIGHT) & 0x1) 
		{
			Player = rotateY(Player, 5000 * delta_time);
		}
		if (GetAsyncKeyState(VK_LEFT) & 0x1) 
		{
			Player = rotateY(Player, -5000 * delta_time);
		}

		if (GetAsyncKeyState('W') & 0x1) 
		{
			//Player[3].xyz += Player[2].xyz * 100 * delta_time;
			defaultView = cameraMove(defaultView, 100, 2, true);
		}
		if (GetAsyncKeyState('S') & 0x1)
		{
			//Player[3].xyz += Player[2].xyz * -100 * delta_time;
			defaultView = cameraMove(defaultView, 100, 2, false);
		}
		if (GetAsyncKeyState('D') & 0x1) 
		{
			//Player = rotateY(Player, 1000 * delta_time);
			defaultView = cameraMove(defaultView, 100, 0, true);
		}
		if (GetAsyncKeyState('A') & 0x1) 
		{
			//Player = rotateY(Player, -1000 * delta_time);
			defaultView = cameraMove(defaultView, 100, 0, false);
		}
		if (GetAsyncKeyState(VK_LBUTTON) & 0x1) 
		{
			//mouseMove(*defaultView);
			//Player = rotateY(Player, -1000 * delta_time);
		}
		//set mouse pos
		if (mbits[0])
		{
			defaultView = mouseMove(defaultView, 1);
		}
		
		prevmPosition = mPosition;

		//Player[3].xyz += Player[0].xyz * 1 * delta_time;
		//Player[3].xyz += Player[2].xyz * 1 * delta_time;
		//Player = rotateY(Player, 100 * delta_time);
		//end::debug_renderer::add_line(Player[3].xyz, Player[3].xyz + Player[0].xyz, Player[0]);
		//end::debug_renderer::add_line(Player[3].xyz, Player[3].xyz + Player[1].xyz, Player[1]);
		//end::debug_renderer::add_line(Player[3].xyz, Player[3].xyz + Player[2].xyz, Player[2]);

		Look_At = (float4x4&)loockATD((XMMATRIX&)Look_At, (XMVECTOR&)Player[3]);
		Turn_To = (float4x4&)turnToD((XMMATRIX&)Turn_To, (XMVECTOR&)Player[3]);

		//Lab 3
		//end::calculate_frustum(frustum, frustumCorners, Player, DirectX::XMConvertToRadians(25), 1.5, 1, 10);
		end::calculate_frustum(frustum, Player, frustumCorners);
		debug_renderer::add_frustum(frustumCorners, B);
		//float4 aabbColor;
		
		for (int i = 0; i < boxesCollision.size(); i++)
		//while (i < boxesCollision.size())
		{
			debug_renderer::add_aabb(boxesCollision[i], 
			(aabb_to_frustum(boxesCollision[i], frustum)) ? G : R);
			//i++;
		}


		//grid
		grid(4);

		//Sorted Emiter 
		pt.pos = { 0, 0, 0 };
		pt.prev_pos = { 0, 0, 0 };
		pt.color = { 1, 1, 1, 1 };
		pt.lifetime = 2;
		//sortedEmiter(pt, 1);

		//Free Emiter 
		fp.spawn_pos = { -5, 0, -5 };
		fp.spawn_color = { 1, 0, 0, 1 };
		//freeEmiter(fp, 1);

		fp2.spawn_pos = { 5, 0, 5 };
		fp2.spawn_color = { 1, 0, 1, 1 };
		//freeEmiter(fp2, 1);

		fp3.spawn_pos = { -5, 0, 5 };
		fp3.spawn_color = { 0, 0, 1, 1 };
		//freeEmiter(fp3, 1);

		fp4.spawn_pos = { 5, 0, -5 };
		fp4.spawn_color = { 0, 1, 0, 1 };
		//freeEmiter(fp4, 1);

		//Matrix Behaviors
		//end::debug_renderer::add_line(shared_pool[pt.indices[i]].prev_pos, shared_pool[pt.indices[i]].pos, shared_pool[pt.indices[i]].color);
		
		return defaultView;
	}
}
