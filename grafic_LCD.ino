#include <TFT_eSPI.h>
#include <Wire.h>
#include <ICM20948_WE.h>
#define ICM20948_ADDR 0x68
#define FOR(i, n) for (int i = 0; i < n; i++)
#define rgb tft.color565

TFT_eSPI tft = TFT_eSPI();
ICM20948_WE myIMU = ICM20948_WE(ICM20948_ADDR);

void applyMatrix(double a, double b, double c, double d,
                 double *x, double *y) {

  // 計算結果を一時的に保持する変数
  double newX;
  double newY;

  // 新しい X 成分 = a*x + b*y
  newX = a * *x + b * *y;

  // 新しい Y 成分 = c*x + d*y
  newY = c * *x + d * *y;

  // 結果を出力ポインタに書き込み
  *x = newX;
  *y = newY;
}

void DrawPixel(int x, int y, uint16_t color) {
  tft.drawPixel(x + tft.width() / 2, -y + tft.height() / 2, color);
}

void DrawLine(int x0, int y0, int x1, int y1, uint16_t color) {
  tft.drawLine(x0 + tft.width() / 2, -y0 + tft.height() / 2, x1 + tft.width() / 2, -y1 + tft.height() / 2, color);
}
void DrawTriangle(int x0, int y0, int x1, int y1, int x2, int y2, uint16_t color) {
  x0 = x0 + tft.width() / 2;
  x1 = x1 + tft.width() / 2;
  x2 = x2 + tft.width() / 2;
  y0 = -y0 + tft.height() / 2;
  y1 = -y1 + tft.height() / 2;
  y2 = -y2 + tft.height() / 2;
  tft.drawTriangle(x0, y0, x1, y1, x2, y2, color);
}
void FillTriangle(int x0, int y0, int x1, int y1, int x2, int y2, uint16_t color) {
  x0 = x0 + tft.width() / 2;
  x1 = x1 + tft.width() / 2;
  x2 = x2 + tft.width() / 2;
  y0 = -y0 + tft.height() / 2;
  y1 = -y1 + tft.height() / 2;
  y2 = -y2 + tft.height() / 2;
  tft.fillTriangle(x0, y0, x1, y1, x2, y2, color);
}

void DrawFixedGUI(int lineSpacing = 5) {

  uint16_t scaleColor = TFT_WHITE;

  DrawLine(-5, 0, 5, 0, scaleColor);  //中心の水平線

  //目盛り生成
  for (int i = 1; i <= 10; i++) {  //上半分
    if (i % 5 == 0) {
      DrawLine(-15, lineSpacing * i, 15, lineSpacing * i, scaleColor);
    } else {
      DrawLine(-10, lineSpacing * i, 10, lineSpacing * i, scaleColor);
    }
  }
  for (int i = 1; i <= 10; i++) {  //下半分
    if (i % 5 == 0) {
      DrawLine(-15, lineSpacing * (-i), 15, lineSpacing * (-i), scaleColor);
    } else {
      DrawLine(-10, lineSpacing * (-i), 10, lineSpacing * (-i), scaleColor);
    }
  }

  //進行方向を示す三角形
  uint16_t dir_tri_color = TFT_WHITE;
  uint16_t dir_tri_outline_clr = TFT_WHITE;
  int ax = 13;
  int bx = ax + 48;
  int cx = bx;
  int dx = cx - 19;
  int ay = 0;
  int by = ay - 15;
  int cy = by - 7;
  int dy = cy;
  DrawTriangle(ax, ay, bx, by, dx, dy, dir_tri_outline_clr);
  DrawTriangle(bx, by, cx, cy, dx, dy, dir_tri_outline_clr);

  DrawTriangle(-ax, ay, -bx, by, -dx, dy, dir_tri_outline_clr);
  DrawTriangle(-bx, by, -cx, cy, -dx, dy, dir_tri_outline_clr);

  //画面の横真っ二つの線がわかるようにある三角形
  uint16_t baseTriangleColor = TFT_WHITE;
  int ex = 100;
  int fx = ex + 25;
  int ey = 0;
  int fy = 10;
  DrawTriangle(ex, ey, fx, fy, fx, -fy, baseTriangleColor);
  DrawTriangle(-ex, ey, -fx, fy, -fx, -fy, baseTriangleColor);


  //ロール目盛りの生成
  uint16_t rollScaleColor = TFT_WHITE;
  const double gy_dis = 200;  //原点からの円弧の距離
  int scaleLength = 10;
  const double scaleIntervalAngle = 5 * (PI / 180);  //度に(pai/180)を掛けるとラジアンに変換できる。
  int tickCount = 7;


  double gx = 0;
  double gy = gy_dis;
  double hx = 0;
  double hy = gy + scaleLength;
  double startAngle = -scaleIntervalAngle * (tickCount + 1);
  applyMatrix(cos(startAngle), -sin(startAngle), sin(startAngle), cos(startAngle), &gx, &gy);
  applyMatrix(cos(startAngle), -sin(startAngle), sin(startAngle), cos(startAngle), &hx, &hy);
  FOR(i, tickCount * 2 + 1) {
    applyMatrix(cos(scaleIntervalAngle), -sin(scaleIntervalAngle), sin(scaleIntervalAngle), cos(scaleIntervalAngle), &gx, &gy);
    applyMatrix(cos(scaleIntervalAngle), -sin(scaleIntervalAngle), sin(scaleIntervalAngle), cos(scaleIntervalAngle), &hx, &hy);
    DrawLine(gx, gy, hx, hy, rollScaleColor);
  }

  //ロール目盛りの円弧
  int dar = gy_dis * sin(scaleIntervalAngle * tickCount);  //一番端っこの目盛りの付け根のx座標
  for (int i = -dar; i < dar; i++) {
    DrawPixel(i, (int)sqrt((gy_dis * gy_dis) - (double)i * i), rollScaleColor);
  }
}


typedef struct {
  float x;
  float y;
  float z;
  float w;
} Vector;

typedef struct {
  Vector a;
  Vector b;
} Line;

#define SIZE 4
void multiply_matrix(const float A[SIZE][SIZE], const float B[SIZE][SIZE], float C[SIZE][SIZE]) {
  // 行列の積 C[i][j] = Sum(A[i][k] * B[k][j])
  for (int i = 0; i < SIZE; i++) {      // 結果Cの行 (i)
    for (int j = 0; j < SIZE; j++) {    // 結果Cの列 (j)
      C[i][j] = 0.0f;                   // 初期化
      for (int k = 0; k < SIZE; k++) {  // 共通のインデックス (k)
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
}

void applyMatrix4(const float M[SIZE][SIZE], Vector *V) {
  Vector temp;
  temp.x = M[0][0] * V->x + M[0][1] * V->y + M[0][2] * V->z + M[0][3] * V->w;
  temp.y = M[1][0] * V->x + M[1][1] * V->y + M[1][2] * V->z + M[1][3] * V->w;
  temp.z = M[2][0] * V->x + M[2][1] * V->y + M[2][2] * V->z + M[2][3] * V->w;
  temp.w = M[3][0] * V->x + M[3][1] * V->y + M[3][2] * V->z + M[3][3] * V->w;
  V->x = temp.x;
  V->y = temp.y;
  V->z = temp.z;
  V->w = temp.w;
}

const double z_max = 1000;
const int h_num = 10;
const double rs = z_max / h_num;  //横線間隔
const int v_num = 10;

Line hline[h_num];
void hlineinitializer() {  //これは直線だからとりあえず長さを1の線分として定義
  for (int i = 0; i < h_num; i++) {
    hline[i].a.x = 0;
    hline[i].b.x = 1;

    hline[i].a.y = hline[i].b.y = 0;

    hline[i].a.z = hline[i].b.z = i * rs;

    hline[i].a.w = hline[i].b.w = 1;
  }
}
Line vline[v_num];  //始点z_max,終点0の線分
void vlineinitializer() {
  for (int i = 0; i < v_num; i++) {
    vline[i].a.x = vline[i].b.x = i * rs;
    vline[i].a.y = vline[i].b.y = 0;
    vline[i].a.z = z_max;
    vline[i].b.z = 0;
    vline[i].a.w = vline[i].b.w = 1;
  }
}


void setup() {
  Wire.begin();
  Serial.begin(115200);
  if (!myIMU.init()) {
    Serial.println("ICM20948 does not respond");
  } else {
    Serial.println("ICM20948 is connected");
  }
  myIMU.autoOffsets();
  myIMU.setAccRange(ICM20948_ACC_RANGE_2G);
  myIMU.setAccDLPF(ICM20948_DLPF_6);
  myIMU.setAccSampleRateDivider(10);

  tft.begin();
  tft.setRotation(0);
  tft.fillScreen(TFT_BLACK);

  //ここで方眼を初期化
  hlineinitializer();
  vlineinitializer();
}

float T[SIZE][SIZE] = {
  { 1, 0, 0, 0 },
  { 0, 1, 0, 0 },
  { 0, 0, 1, 0 },
  { 0, 0, 0, 1 }
} 
float Rx[SIZE][SIZE] = {
  { 1, 0, 0, 0 },
  { 0, 1, 0, 0 },
  { 0, 0, 1, 0 },
  { 0, 0, 0, 1 }
} 
float Rz[SIZE][SIZE] = {
  { 1, 0, 0, 0 },
  { 0, 1, 0, 0 },
  { 0, 0, 1, 0 },
  { 0, 0, 0, 1 }
} 
float M_temp[SIZE][SIZE];
float M[SIZE][SIZE];

int y0a = 0;
int y1a = 0;
void loop() {

  myIMU.readSensor();

  float pitch = myIMU.getPitch();
  float roll = myIMU.getRoll();
  Serial.print("Pitch:  ");
  Serial.print(pitch, 10);
  Serial.print("  |  Roll: ");
  Serial.print(roll, 10);
  Serial.println("");

  float Spacing = 10;

  float h = 10;  //高度単位は画面px
  T[1][3] = -h;  //下にh文下げる

  float fai = -pitch * PI / 180;
  Rx[1][1] = cos(fai);Rx[1][2] = sin(fai);
  Rx[2][1] = -sin(fai);Rx[2][2] = cos(fai)

  float theta = -roll * PI / 180.0;
  Rz[0][0] = cos(theta);Rz[1][0] = sin(theta);
  Rz[0][1] = -sin(theta);Rz[1][1] = cos(theta);

  multiply_matrix(Rz,Rx,M_temp);
  multiply_matrix(M_temp,T,M);

  for (int i = 0; i < h_num; i++) {
    applyMatrix4(M,hline[i].a);
    applyMatrix4(M,hline[i].b);
  }
  for (int i = 0; i < v_num; i++) {
    applyMatrix4(M,vline[i].a);
    applyMatrix4(M,vline[i].b);
  }
  DrawLine(-160, y0a, 160, y1a, TFT_BLACK);
  y0a = tan(theta) * (-160) - Spacing * pitch;
  y1a = tan(theta) * (160) - Spacing * pitch;
  DrawLine(-160, y0a, 160, y1a, TFT_GREEN);

  DrawFixedGUI(Spacing);
}