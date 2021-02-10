#include <Triangle.h>

float Triangle::perim() {
  return this->node1.distance(this->node2) + \
         this->node2.distance(this->node3) + \
         this->node3.distance(this->node1);
}
float Triangle::longest() {
  float a = this->node1.distance(this->node2), b = this->node2.distance(this->node3), c = this->node3.distance(this->node1), longest;
  float temp = (a >= b ? a : b);
  return (c >= temp ? c : temp);
}
float Triangle::shortest() {
  float a = this->node1.distance(this->node2), b = this->node2.distance(this->node3), c = this->node3.distance(this->node1), longest;
  float temp = (a <= b ? a : b);
  return (c <= temp ? c : temp);
}
float Triangle::area() {
  float p = this->perim()/2;
  return sqrt(p* \
    (p-this->node1.distance(this->node2))* \
    (p-this->node2.distance(this->node3))* \
    (p-this->node3.distance(this->node1)));
}
float Triangle::cx() {
  return (this->node1.x+this->node2.x+this->node3.x)/3;
}
float Triangle::cy() {
  return (this->node1.y+this->node2.y+this->node3.y)/3;
}
float Triangle::cond() {
  return this->area()/(this->perim());
}
