#include <Node.h>

double Node::distance(Node obj) {
  return sqrt((this->x-obj.x)*(this->x-obj.x)+(this->y-obj.y)*(this->y-obj.y));
}
