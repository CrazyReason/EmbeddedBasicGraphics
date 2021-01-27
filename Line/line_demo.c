#include "svpng.inc"
#include <math.h> // fminf(), fmaxf(), sinf(), cosf()
#define PI 3.14159265359f
#define W 512
#define H 512
static unsigned char img[W * H * 3];

#define ipart_(X) ((int)(X))
#define round_(X) ((int)(((float)(X))+0.5f))
#define fpart_(X) (((float)(X))-(float)ipart_(X))
#define rfpart_(X) (1.0-fpart_(X))

#define swap_(a,b)    (a=(a)+(b),b=(a)-(b),a=(a)-(b))



void setPixel(int x, int y) {
	unsigned char* p = img + (y * W + x) * 3;
	p[0] = p[1] = p[2] = 255;
}

void alphablend(int x, int y, float alpha, float r, float g, float b) {
	unsigned char* p = img + (y * W + x) * 3;
	p[0] = (unsigned char)(p[0] * (1 - alpha) + r * alpha );
	p[1] = (unsigned char)(p[1] * (1 - alpha) + g * alpha );
	p[2] = (unsigned char)(p[2] * (1 - alpha) + b * alpha );
}

void setPixelAlpha(int x, int y, float brightness)
{
	alphablend(x, y, brightness, 255.0f, 255.0f, 255.0f);
}


//缺点 还是不算快 存在浮点计算 而且大量浮点累加 可能存在误差
void lineDDA(int x0, int y0, int x1, int y1)
{
	int dx = x1 - x0, dy = y1 - y0, steps, k;
	float xIncrement, yIncrement, x = x0, y = y0;

	if (fabs(dx) > fabs(dy))	//判断增长方向
		steps = fabs(dx);		//以X为单位间隔取样计算
	else
		steps = fabs(dy);		//以Y为单位间隔取样计算

	xIncrement = (float)(dx) / (float)(steps);	//计算每个度量间隔的X方向增长量
	yIncrement = (float)(dy) / (float)(steps);	//计算每个度量间隔的Y方向增长量

	while(steps--)
	{
		setPixel(round_(x), round_(y));
		x += xIncrement;
		y += yIncrement;
		
	}
}

static inline setPixelHorWidth(int x, int y, int w)
{
	while (w--)
	{
		setPixel(x + w, y);
	}
}

static inline setPixeVerWidth(int x, int y, int w)
{
	while (w--)
	{
		setPixel(x, y + w);
	}
}

//Bresenham 算法
void lineBres(int x0, int y0, int x1, int y1,int r)
{
	int p, twoDy, twoDyMinusDx, s1, s2;
	int dx = abs(x1 - x0), dy = abs(y1 - y0);

	if (dy > dx)	//斜率大于1 
	{
		p = 2 * dx - dy;
		twoDy = 2 * dx;
		twoDyMinusDx = 2 * (dx - dy);

		if (y0 > y1)//斜率为负时 反转斜率
		{
			swap_(x0, x1);
			swap_(y0, y1);
		}
		s1 = x1 > x0 ? 1 : -1;

		setPixelHorWidth(x0, y0,r);

		while (y0 < y1)
		{
			y0++;
			if (p < 0)
			{
				p += twoDy;
			}
			else
			{
				x0 += s1;
				p += twoDyMinusDx;
			}
			setPixelHorWidth(x0, y0, r);
		}
	}
	else
	{
		p = 2 * dy - dx;
		twoDy = 2 * dy;
		twoDyMinusDx = 2 * (dy - dx);

		if (x0 > x1)//斜率为负时 反转斜率
		{
			swap_(x0, x1);
			swap_(y0, y1);
		}
		s2 = y1 > y0 ? 1 : -1;

		setPixeVerWidth(x0, y0, r);

		while (x0 < x1)
		{
			x0++;
			if (p < 0)
			{
				p += twoDy;
			}
			else
			{
				y0 += s2;
				p += twoDyMinusDx;
			}
			setPixeVerWidth(x0, y0, r);
		}
	}
}


//抗锯齿算法
void lineAnti_Wu(int x0, int y0, int x1, int y1)
{
	int steep = abs(y1 - y0) > abs(x1 - x0);

	// swap the co-ordinates if slope > 1 or we 
	// draw backwards 
	if (steep)
	{
		swap_(x0, y0);
		swap_(x1, y1);
	}
	if (x0 > x1)
	{
		swap_(x0, x1);
		swap_(y0, y1);
	}

	//compute the slope 
	float dx = x1 - x0;
	float dy = y1 - y0;
	float gradient = dy / dx;
	if (dx == 0.0)
		gradient = 1;

	int xpxl1 = x0;
	int xpxl2 = x1;
	float intersectY = y0;
	int  x;

	// main loop 
	if (steep)
	{
		for (x = xpxl1; x <= xpxl2; x++)
		{
			setPixelAlpha(ipart_(intersectY), x, rfpart_(intersectY));
			setPixelAlpha(ipart_(intersectY) + 1, x,fpart_(intersectY));
			intersectY += gradient;
		}
	}
	else
	{
		for (x = xpxl1; x <= xpxl2; x++)
		{
			setPixelAlpha(x, ipart_(intersectY),rfpart_(intersectY));
			setPixelAlpha(x, ipart_(intersectY)+1,fpart_(intersectY));
			intersectY += gradient;
		}
	}
}

void lineAnti_WuMulti(int x0, int y0, int x1, int y1, int r)
{
	int steep = abs(y1 - y0) > abs(x1 - x0);

	// swap the co-ordinates if slope > 1 or we 
	// draw backwards 
	if (steep)
	{
		swap_(x0, y0);
		swap_(x1, y1);
	}
	if (x0 > x1)
	{
		swap_(x0, x1);
		swap_(y0, y1);
	}

	//compute the slope 
	float dx = x1 - x0;
	float dy = y1 - y0;
	float gradient = dy / dx;
	if (dx == 0.0)
		gradient = 1;

	int xpxl1 = x0;
	int xpxl2 = x1;
	float intersectY = y0;

	int i, x;

	// main loop 
	if (steep)
	{
		for (x = xpxl1; x <= xpxl2; x++)
		{
			// pixel coverage is determined by fractional 
			// part of y co-ordinate 		
			setPixelAlpha(ipart_(intersectY), x, rfpart_(intersectY));
			for (i = 1; i < r; i++)
			{
				setPixel(ipart_(intersectY) + i, x);
			}
			setPixelAlpha(ipart_(intersectY) + r, x, fpart_(intersectY));
			intersectY += gradient;
		}
	}
	else
	{
		for (x = xpxl1; x <= xpxl2; x++)
		{
			// pixel coverage is determined by fractional 
			// part of y co-ordinate 		
			setPixelAlpha(x, ipart_(intersectY), rfpart_(intersectY));
			for (i = 1; i < r; i++)
			{
				setPixel(x, ipart_(intersectY) + i);
			}
			setPixelAlpha(x, ipart_(intersectY) + r, fpart_(intersectY));
			intersectY += gradient;
		}
	}
}

void lineAnti_AreaWeight(int x0, int y0, int x1, int y1, int r)
{
	//高斯核进行权重划分  边缘锐化效果
#if 1
	//int weight[5][5] = { { 1, 2, 4, 2, 1 }, { 2, 5,6,5, 2 }, { 4, 6, 8,6,4 },{2,5,6,5,2 },{1,2,4,2,1} };
	//使用高斯累计核，优化速度
	int weight_sum[5][6] = { {0, 1, 3, 7, 9, 10 }, {0, 2, 7,13,18, 20 }, {0, 4, 10, 18,24,28 },{0,2, 7,13,18, 20 },{0,1, 3, 7, 9, 10} };
#else
	//使用平均核 查看效果  边缘平滑
	//int weight[5][5] = { { 3, 3, 5, 3, 3 }, { 3, 3, 5, 3, 3 }, { 4, 4, 4, 4, 4 },{ 3, 3, 5, 3, 3 },{ 3, 3, 5, 3, 3 } };
	int weight_sum[5][6] = { {0, 3, 6, 11, 14, 17 }, {0, 3, 6, 11, 14, 17 }, {0, 4, 8, 12,16,20 },{0, 3, 6, 11, 14, 17 },{0, 3, 6, 11, 14, 17 } };
#endif

	int weight_1, weight_2, weight_3;
	float weight_temp;

	int steep = abs(y1 - y0) > abs(x1 - x0);

	// swap the co-ordinates if slope > 1 or we 
	// draw backwards 
	if (steep)
	{
		swap_(x0, y0);
		swap_(x1, y1);
	}
	if (x0 > x1)
	{
		swap_(x0, x1);
		swap_(y0, y1);
	}

	//compute the slope 
	float dx = x1 - x0;
	float dy = y1 - y0;
	float gradient = dy / dx;
	if (dx == 0.0f)
		gradient = 1;

	int xpxl1 = x0;
	int xpxl2 = x1;
	float intersectY = y0;

	int i, j, x;
	float temp;

	//上边界方程 y-y0 = k*(x -x0)
	//下边界方程 y-y0 = k*(x -x0)-r

	//因0<K<1 每条边界线只存在两种情况需要透明度计算
	//1.穿过一个像素
	//2.穿过两个像素
	//由此可知R+2为最大像素宽度，不需要抗锯齿的像素位置为 (Y-R +1) < int(Y) < (Y + 1 -2)
	//则边界线之间的距离小于3，则还会存在同时穿过一个像素

	// main loop 
	if (steep)
	{
		for (x = xpxl1; x <= xpxl2; x++)
		{
			weight_1 = 0;
			weight_2 = 0;
			weight_3 = 0;

			weight_temp = intersectY;
			for (j = 0; j < 5; j++)
			{
				//上边界经过的像素点
				weight_temp += gradient / 5;
				temp = fabs(weight_temp - ipart_(intersectY)) * 5;
				if (temp > 5.0f)
				{
					weight_2 += weight_sum[j][ipart_(temp - 5.0f)];
					weight_1 += weight_sum[j][5];	//当穿越时，下面的像素应该是满强度的

					//下边界经过的像素点 由直线的对称性得出
					weight_3 += weight_sum[j][0];
				}
				else
				{
					weight_1 += weight_sum[j][ipart_(temp)];

					//下边界经过的像素点 由直线的对称性得出
					weight_3 += weight_sum[j][ipart_(5.0f - temp)];
				}
			}
			// pixel coverage is determined by fractional 
			// part of y co-ordinate 
			if (temp > 5.0f)
			{
				setPixelAlpha(ipart_(intersectY) + 1, x, (88 - weight_2) / 88.0f);
				setPixelAlpha(ipart_(intersectY), x, weight_3 / 88.0f);

				setPixelAlpha(ipart_(intersectY) + r, x, weight_1 / 88.0f);
				setPixelAlpha(ipart_(intersectY) + r + 1, x, weight_2 / 88.0f);

				for (i = 1; i < r; i++)
				{
					setPixel(ipart_(intersectY) + i, x);
				}
			}
			else
			{
				setPixelAlpha(ipart_(intersectY) + r, x, weight_1 / 88.0f);
				setPixelAlpha((ipart_(intersectY)), x, weight_3 / 88.0f);

				for (i = 1; i < r; i++)
				{
					setPixel(ipart_(intersectY) + i, x);
				}
			}
			intersectY += gradient;
		}
	}
	else
	{
		for (x = xpxl1; x <= xpxl2; x++)
		{
			weight_1 = 0;
			weight_2 = 0;
			weight_3 = 0;
			weight_temp = intersectY;
			for (j = 0; j < 5; j++)
			{
				//上边界经过的像素点
				weight_temp += gradient / 5;
				temp = fabs((weight_temp - ipart_(intersectY)) * 5);
				if (temp > 5.0f)
				{
					weight_2 += weight_sum[j][ipart_(temp - 5.0f)];
					weight_1 += weight_sum[j][5];	//当穿越时，下面的像素应该是满强度的

					//下边界经过的像素点 由直线的对称性得出
					weight_3 += weight_sum[j][0];
				}
				else
				{
					weight_1 += weight_sum[j][ipart_(temp)];

					//下边界经过的像素点 由直线的对称性得出
					weight_3 += weight_sum[j][ipart_(5.0f - temp)];
				}
			}

			// pixel coverage is determined by fractional 
			// part of y co-ordinate 
			if (temp > 5.0f)
			{
				setPixelAlpha(x, ipart_(intersectY), weight_3 / 88.0f);
				setPixelAlpha(x, ipart_(intersectY) + 1, (88 - weight_2) / 88.0f);

				setPixelAlpha(x, ipart_(intersectY) + r, weight_1 / 88.0f);
				setPixelAlpha(x, ipart_(intersectY) + r + 1, weight_2 / 88.0f);

				for (i = 2; i < r; i++)
				{
					setPixel(x, ipart_(intersectY) + i);
				}
			}
			else
			{
				setPixelAlpha(x, ipart_(intersectY), weight_3 / 88.0f);
				setPixelAlpha(x, (ipart_(intersectY)) + r, weight_1 / 88.0f);

				for (i = 1; i < r; i++)
				{
					setPixel(x, ipart_(intersectY) + i);
				}
			}
			intersectY += gradient;
		}
	}
}


//#define Line(X0,Y0,X1,Y1,R)		lineDDA(X0,Y0,X1,Y1)
//#define Line(X0,Y0,X1,Y1,R)		lineBres(X0,Y0,X1,Y1,R)
//#define Line(X0,Y0,X1,Y1,R)			lineAnti_Wu(X0,Y0,X1,Y1)
#define Line(X0,Y0,X1,Y1,R)			lineAnti_WuMulti(X0,Y0,X1,Y1,R)
//#define Line(X0,Y0,X1,Y1,R)			lineAnti_AreaWeight(X0,Y0,X1,Y1,R)


int main() {

	memset(img, 0, sizeof(img));

	float cx = W * 0.5f, cy = H * 0.5f;
	for (int j = 0; j < 5; j++) {
		float r1 = fminf(W, H) * (j + 0.5f) * 0.085f;
		float r2 = fminf(W, H) * (j + 1.5f) * 0.085f;
		float t = j * PI / 64.0f, r = (j + 1) * 1.0f;
		for (int i = 1; i <= 64; i++, t += 2.0f * PI / 64.0f) {
			float ct = cosf(t), st = sinf(t);
			Line(cx + r1 * ct, cy - r1 * st, cx + r2 * ct, cy - r2 * st, r);
		}
	}

	//Line(100, 100, 200, 150, 10);

	svpng(fopen("line_WuMulti.png", "wb"), W, H, img, 0);
}


