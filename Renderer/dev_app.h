#pragma once
#include <cstdint>
#include <chrono>
#include "renderer.h"

namespace end
{
	// Simple app class for development and testing purposes
	struct dev_app_t
	{
		//void update();
		view_t update(view_t defaultView);
		dev_app_t();

		void setMDown(uint8_t number);
		void setMUp(uint8_t number);
		void setMPosition(float x, float y);

		double get_delta_time()const;
	};
}