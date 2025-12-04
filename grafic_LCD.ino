#include <TFT_eSPI.h>
TFT_eSPI tft = TFT_eSPI(); 
void setup() {
  Serial.begin(115200);
  Serial.println("wtf");
  tft.init();
  tft.setRotation(1); 
  tft.fillScreen(TFT_BLACK); 
  tft.setTextColor(TFT_WHITE, TFT_BLACK); 
  tft.drawString("Hello ILI9488 World!", 50, 100, 4); 
}
void loop() {
}