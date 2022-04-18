#include "debug_renderer.h"
#include <array>

// Anonymous namespace
namespace
{
	// Declarations in an anonymous namespace are global BUT only have internal linkage.
	// In other words, these variables are global but are only visible in this source file.

	// Maximum number of debug lines at one time (i.e: Capacity)
	constexpr size_t MAX_LINE_VERTS = 4096*4; 

	// CPU-side buffer of debug-line verts
	// Copied to the GPU and reset every frame.
	size_t line_vert_count = 0;
	std::array< end::colored_vertex, MAX_LINE_VERTS> line_verts;

	//Matrices
	
}

namespace end
{
	namespace debug_renderer
	{
		float3 point;
		float3 point2;
		float3 point3;
		float3 point4;
		float3 point5;
		float3 point6;
		float3 point7;
		float3 point8;

		//TODO Once you finish thisfile correctly, you should see a green checkmark when you run the game.
		void add_line(float3 point_a, float3 point_b, float4 color_a, float4 color_b)
		{
			//TODO Add points to debug_verts, increments debug_vert_count
			//line_verts = {  point_a,  point_b,  color_a,  color_b };
			//get_line_verts()->pos = { point_a };
			//const colored_vertex a;
			//a.pos = { point_a };
			/*for (size_t i = 0; i < line_verts.size(); i++)
			{
				line_verts[i].pos = point_a;
				line_verts[i].pos = point_b;
				line_verts[i].color = color_a;
				line_verts[i].color = color_b;
			}*/
			line_verts[line_vert_count].pos = point_a;
			line_verts[line_vert_count].color = color_a;
			line_vert_count++;
			line_verts[line_vert_count].pos = point_b;
			line_verts[line_vert_count].color = color_b;
			line_vert_count++;

		}

		void clear_lines()
		{
			//TODO Resets debug_vert_count to 0
			line_vert_count = 0;
		}

		const colored_vertex* get_line_verts()
		{ 
			//TODO Returns the line vert array pointer
			std::array< end::colored_vertex, MAX_LINE_VERTS> LV = line_verts;
			return LV._Elems;
		}

		size_t get_line_vert_count() 
		{ 
			//TODO Returns how many vertices there are now

			return line_vert_count;
		}

		size_t get_line_vert_capacity()
		{
			//TODO returns the maximum size the line vert array

			return MAX_LINE_VERTS;
		}

		float3 centerEqu(aabb_t aabb, bool check)
		{
			float3 point;
			if (check)
			{
				point = aabb.center + aabb.extents;
			}
			else
			{
				point = aabb.center - aabb.extents;
			}

			return point;
		}
		float3 extentsEqu(aabb_t aabb)
		{
			float3 point;
			point.x -= 2 * aabb.extents.x;
			return point;
		}
		void add_aabbExtend(aabb_t aabb)
		{
			//
			point3.x -= 2 * aabb.extents.x;
			point2.y -= 2 * aabb.extents.y;
			point4.z -= 2 * aabb.extents.z;
			//
			point5.y += 2 * aabb.extents.y;
			point6.x += 2 * aabb.extents.x;
			point7.z += 2 * aabb.extents.z;
		}

		void add_aabbCenter(aabb_t aabb)
		{
			//front
			point = centerEqu(aabb, true);
			point2 = centerEqu(aabb, true);
			point3 = centerEqu(aabb, true);
			point4 = centerEqu(aabb, true);
			//back
			point5 = centerEqu(aabb, false);
			point6 = centerEqu(aabb, false);
			point7 = centerEqu(aabb, false);
			point8 = centerEqu(aabb, false);
		}
		
		void add_aabb(aabb_t aabb, float4 color)
		{
			
			add_aabbCenter(aabb);
			add_aabbExtend( aabb);
			
			//
			add_line(point, point2, color, color);
			add_line(point, point3, color, color);
			add_line(point, point4, color, color);
			add_line(point8, point5, color, color);
			add_line(point8, point6, color, color);
			add_line(point8, point7, color, color);
			//
			add_line(point4, point6, color, color);
			add_line(point3, point7, color, color);
			add_line(point2, point6, color, color);
			add_line(point2, point7, color, color);
			add_line(point5, point3, color, color);
			add_line(point5, point4, color, color);
		}

		void add_frustum(points f_vertices, float4 color)
		{
			// far plane
			//for (size_t i = 0; i < 3; i++)
			//{
			//	//add_line(f_vertices[i], f_vertices[i], color, color);
			//}
			add_line(f_vertices[0], f_vertices[1], color, color);
			add_line(f_vertices[2], f_vertices[0], color, color);
			add_line(f_vertices[3], f_vertices[1], color, color);
			
			add_line(f_vertices[3], f_vertices[2], color, color);

			// conecting lines
			for (size_t i = 0; i < 4; i++)
			{
				add_line(f_vertices[i], f_vertices[i + 4], color, color);
			}
		
			// near plane
			add_line(f_vertices[4], f_vertices[5], color, color);
			add_line(f_vertices[4], f_vertices[6], color, color);
			add_line(f_vertices[5], f_vertices[7], color, color);
			add_line(f_vertices[6], f_vertices[7], color, color);


		}
	}

}