#include "image.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>

// allows the use of min() and max() functions
#include <algorithm>
using namespace std;


Pixel::Pixel(const Pixel32 &p) {
	r = p.r / 255.0;
	g = p.g / 255.0;
	b = p.b / 255.0;
	a = p.a / 255.0;
}

Pixel32::Pixel32(const Pixel &p) {
	r = p.r <= 1 ? (p.r >= 0 ? p.r * 255 : 0) : 255;
	g = p.g <= 1 ? (p.g >= 0 ? p.g * 255 : 0) : 255;
	b = p.b <= 1 ? (p.b >= 0 ? p.b * 255 : 0) : 255;
	a = p.a <= 1 ? (p.a >= 0 ? p.a * 255 : 0) : 255;
}

bool inBounds(int x, int y, int w, int h) {
	return x >= 0 && x < w && y >= 0 && y < h;
}

int Image32::AddRandomNoise(const float &noise, Image32 &outputImage) const {
	srand(time(NULL));
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			p.r += (((rand() % 1000) - 500) / 500.0) * noise;
			p.g += (((rand() % 1000) - 500) / 500.0) * noise;
			p.b += (((rand() % 1000) - 500) / 500.0) * noise;
			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}

int Image32::Brighten(const float &brightness, Image32 &outputImage) const {
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			p.r *= brightness;
			p.g *= brightness;
			p.b *= brightness;
			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}

int Image32::Luminance(Image32 &outputImage) const {
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			float l = 0.3 * p.r + 0.59 * p.g + 0.11 * p.b;
			p.r = l;
			p.g = l;
			p.b = l;
			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}

int Image32::Contrast(const float &contrast, Image32 &outputImage) const {
	float luminance = 0.0;
	int pix = w * h;
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			luminance += 0.3 * p.r + 0.59 * p.g + 0.11 * p.b;
		}
	}
	luminance = luminance / pix;
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			p.r = luminance + (p.r - luminance) * contrast;
			p.g = luminance + (p.g - luminance) * contrast;
			p.b = luminance + (p.b - luminance) * contrast;
			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}

int Image32::Saturate(const float &saturation, Image32 &outputImage) const {
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			float luminance = 0.3 * p.r + 0.59 * p.g + 0.11 * p.b;
			p.r = luminance + (p.r - luminance) * saturation;
			p.g = luminance + (p.g - luminance) * saturation;
			p.b = luminance + (p.b - luminance) * saturation;
			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}

int Image32::Quantize(const int &bits, Image32 &outputImage) const {
	int maximum = 0xff >> (8 - bits);
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			p.r = (float) ((int) (p.r * (maximum + 1))) / maximum;
			p.g = (float) ((int) (p.g * (maximum + 1))) / maximum;
			p.b = (float) ((int) (p.b * (maximum + 1))) / maximum;
			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}

int Image32::RandomDither(const int &bits, Image32 &outputImage) const {
	AddRandomNoise(0.3, outputImage);
	Quantize(bits, outputImage);
	return 1;
}

int Image32::OrderedDither2X2(const int &bits, Image32 &outputImage) const {
	int maximum = 0xff >> (8 - bits);
	int dith[4] = {1, 3, 4, 2};
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			int index = (i % 2) + (j % 2) * 2;
			Pixel p = Pixel(outputImage.pixel(i, j));
			p.r = (float) ((int) (p.r * maximum + (1 - dith[index] / 5.0))) / maximum;
			p.g = (float) ((int) (p.g * maximum + (1 - dith[index] / 5.0))) / maximum;
			p.b = (float) ((int) (p.b * maximum + (1 - dith[index] / 5.0))) / maximum;
			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}

int Image32::FloydSteinbergDither(const int &bits, Image32 &outputImage) const {
	int maximum = 0xff >> (8 - bits);
	vector< vector<float> > rErrors(w + 1, vector<float>(h + 1, 0));
	vector< vector<float> > gErrors(w + 1, vector<float>(h + 1, 0));
	vector< vector<float> > bErrors(w + 1, vector<float>(h + 1, 0));
	for (int j = 0; j < h; j++) {
		for (int i = 0; i < w; i++) {
			Pixel p = Pixel(outputImage.pixel(i, j));
			float rError = rErrors[i][j];
			float gError = gErrors[i][j];
			float bError = bErrors[i][j];

			p.r += rError;
			float r = (float) ((int) (p.r * (maximum + 1))) / maximum;
			float rE = p.r - r;
			p.r = r;

			p.g += gError;
			float g = (float) ((int) (p.g * (maximum + 1))) / maximum;
			float gE = p.g - g;
			p.g = g;

			p.b += bError;
			float b = (float) ((int) (p.b * (maximum + 1))) / maximum;
			float bE = p.b - b;
			p.b = b;

			outputImage.pixel(i, j) = Pixel32(p);

			rErrors[i + 1][j] += rE * 0.4375;
			rErrors[i][j + 1] += rE * 0.3125;
			rErrors[i + 1][j + 1] += rE * 0.0625;
			if (i != 0) rErrors[i - 1][j + 1] += rE * 0.1875;

			gErrors[i + 1][j] += gE * 0.4375;
			gErrors[i][j + 1] += gE * 0.3125;
			gErrors[i + 1][j + 1] += gE * 0.0625;
			if (i != 0) gErrors[i - 1][j + 1] += gE * 0.1875;

			bErrors[i + 1][j] += bE * 0.4375;
			bErrors[i][j + 1] += bE * 0.3125;
			bErrors[i + 1][j + 1] += bE * 0.0625;
			if (i != 0) bErrors[i - 1][j + 1] += bE * 0.1875;
		}
	}
	return 1;
}

int Image32::Blur3X3(Image32 &outputImage) const {
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			outputImage.pixel(i, j) = GaussianSample(i, j, 0.7, 1.5);
		}
	}
	return 1;
}

int Image32::EdgeDetect3X3(Image32 &outputImage) const {
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			Pixel p1 = Pixel();

			int weight = 0;
			for (int p = -1; p <= 1; ++p) {
				for (int q = -1; q <= 1; ++q) {
					int x = i + p;
					int y = j + q;
					if (inBounds(x, y, w, h)) {
						Pixel p2 = Pixel(pixel(x, y));
						if (p != 0 && q != 0) {
							p1.r += p2.r * -1;
							p1.g += p2.g * -1;
							p1.b += p2.b * -1;
							weight++;
						}
					}
				}
			}
			Pixel p2 = Pixel(pixel(i, j));
			p1.r += p2.r * weight;
			p1.g += p2.g * weight;
			p1.b += p2.b * weight;
			outputImage.pixel(i, j) = Pixel32(p1);
		}
	}
	return 1;
}
int Image32::ScaleNearest(const float &scaleFactor, Image32 &outputImage) const {
	int width, height;
	width = w * scaleFactor;
	height = h * scaleFactor;
	outputImage.setSize(width, height);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			outputImage.pixel(i, j) = NearestSample(i / scaleFactor, j / scaleFactor);
		}
	}
	return 1;
}

int Image32::ScaleBilinear(const float &scaleFactor, Image32 &outputImage) const {
	int width, height;
	width = w * scaleFactor;
	height = h * scaleFactor;
	outputImage.setSize(width, height);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			outputImage.pixel(i, j) = BilinearSample(i / scaleFactor, j / scaleFactor);
		}
	}
	return 1;
}

int Image32::ScaleGaussian(const float &scaleFactor, Image32 &outputImage) const {
	int width, height;
	width = w * scaleFactor;
	height = h * scaleFactor;
	outputImage.setSize(width, height);
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			outputImage.pixel(i, j) = GaussianSample(i / scaleFactor, j / scaleFactor, 0.7, 2.5);
		}
	}
	return 1;
}

int Image32::RotateNearest(const float &angle, Image32 &outputImage) const {
	return 0;
}

int Image32::RotateBilinear(const float &angle, Image32 &outputImage) const {
	return 0;
}

int Image32::RotateGaussian(const float &angle, Image32 &outputImage) const {
	return 0;
}

int Image32::SetAlpha(const Image32 &matte) {
	return 0;
}

int Image32::Composite(const Image32 &overlay, Image32 &outputImage) const {
	return 0;
}

int Image32::CrossDissolve(const Image32 &source, const Image32 &destination, const float &blendWeight, Image32 &outputImage) {
	for (int i = 0; i < source.width(); i++) {
		for (int j = 0; j < source.height(); j++) {
			Pixel p1 = Pixel(source.pixel(i, j));
			Pixel p2 = Pixel(destination.pixel(i, j));
			Pixel p = Pixel();

			p.r = blendWeight * p1.r + (1 - blendWeight) * p2.r;
			p.g = blendWeight * p1.g + (1 - blendWeight) * p2.g;
			p.b = blendWeight * p1.b + (1 - blendWeight) * p2.b;

			outputImage.pixel(i, j) = Pixel32(p);
		}
	}
	return 1;
}
int Image32::Warp(const OrientedLineSegmentPairs &olsp, Image32 &outputImage) const {
	outputImage.setSize(w, h);
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++) {
			float x, y;
			olsp.getSourcePosition(i, j, x, y);
			if (x < 0) x = 0;
			else if (x >= w) x = w - 1;
			if (y < 0) y = 0;
			else if (y >= h) y = h - 1;

			outputImage.pixel(i, j) = BilinearSample(x, y);
		}
	}
	return 1;
}

int Image32::FunFilter(Image32 &outputImage) const {
	return 0;
}

int Image32::Crop(const int &x1, const int &y1, const int &x2, const int &y2, Image32 &outputImage) const {
	outputImage.setSize(x2 - x1, y2 - y1);
	for (int i = x1; i < x2; ++i) {
		for (int j = y1; j < y2; ++j) {
			outputImage.pixel(i - x1, j - y1) = this->pixel(i, j);
		}
	}
	return 1;
}

Pixel32 Image32::NearestSample(const float &x, const float &y) const {
	int i = x;
	int j = y;
	if (x - i > 0.5) {
		i += 1;
	}
	if (y - j > 0.5) {
		j += 1;
	}
	return pixel(i, j);
}

Pixel32 Image32::BilinearSample(const float &x, const float &y) const {
	int i = (int) x;
	int j = (int) y;

	float a = x - i;
	float b = y - j;

	float q1 = a * b;
	float q2 = (1 - a) * b;
	float q3 = a * (1 - b);
	float q4 = (1 - a) * (1 - b);

	Pixel p1 = Pixel(pixel(i, j));
	Pixel p2 = i + 1 >= w ? p1 : Pixel(pixel(i + 1, j));
	Pixel p3 = j + 1 >= h ? p1 : Pixel(pixel(i, j + 1));
	Pixel p4 = i + 1 >= w || j + 1 >= h ? p1 : Pixel(pixel(i + 1, j + 1));

	Pixel p = Pixel();

	p.r = q4 * p1.r + q3 * p2.r + q2 * p3.r + q1 * p4.r;
	p.g = q4 * p1.g + q3 * p2.g + q2 * p3.g + q1 * p4.g;
	p.b = q4 * p1.b + q3 * p2.b + q2 * p3.b + q1 * p4.b;

	return Pixel32(p);
}

static const float inv_sqrt_2pi = 0.3989422804014327;
float evalNormal(float d, float var) {
	float a = d / var;
	return inv_sqrt_2pi / sqrt(var) * exp(-0.5 * a * a);
}

float dist(float x1, float y1, float x2, float y2) {
	float dx = x2 - x1;
	float dy = y2 - y1;
	return sqrt(dx * dx + dy * dy);
}


Pixel32 Image32::GaussianSample(const float &x, const float &y, const float &variance, const float &radius) const {
	int x1 = x;
	int y1 = y;
	int rad = radius;
	int diameter = (rad + 1) * 2;
	float totalWeight = 0;
	Pixel p = Pixel();
	for (int i = 0; i < diameter; ++i) {
		for (int j = 0; j < diameter; ++j) {
			int a = x1 + i - rad;
			int b = y1 + j - rad;
			if (inBounds(a, b, w, h)) {
				float distance = dist((float) a, (float) b, x, y);
				if (distance <= radius) {
					Pixel p2 = Pixel(pixel(a, b));
					float weight = evalNormal(distance, variance);
					totalWeight += weight;
					p.r += p2.r * weight;
					p.g += p2.g * weight;
					p.b += p2.b * weight;
				}
			}
		}
	}
	p.r = p.r / totalWeight;
	p.g = p.g / totalWeight;
	p.b = p.b / totalWeight;
	return Pixel32(p);
}
