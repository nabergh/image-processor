#include "lineSegments.h"
#include <math.h>

float pyth(float x1, float y1, float x2, float y2) {
	float dx = x2 - x1;
	float dy = y2 - y1;
	return sqrt(dx * dx + dy * dy);
}

float OrientedLineSegment::length(void) const {
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

float OrientedLineSegment::distance(const int &x, const int &y) const {
	float sX = x - x1;
	float sY = y - y1;
	float pX, pY;
	getPerpendicular(pX, pY);
	return sX * pX + sY * pY;
}

void  OrientedLineSegment::getPerpendicular(float &x, float &y) const {
	float l = length();
	x = (y2 - y1) / l;
	y = -1 * (x2 - x1) / l;
}

void  OrientedLineSegment::GetSourcePosition(const OrientedLineSegment &source, const OrientedLineSegment &destination,
        const int &targetX, const int &targetY,
        float &sourceX, float &sourceY) {
	float length = destination.length();
	float top = (targetX - destination.x1) * (destination.x2 - destination.x1) + (targetY - destination.y1) * (destination.y2 - destination.y1);
	float u = top / (length * length);
	float vX = source.x1 + u * (source.x2 - source.x1);
	float vY = source.y1 + u * (source.y2 - source.y1);

	float pX, pY;
	source.getPerpendicular(pX, pY);
	float dist = destination.distance(targetX, targetY);
	pX = pX * dist;
	pY = pY * dist;

	sourceX = vX + pX;
	sourceY = vY + pY;
}
