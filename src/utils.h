#ifndef UTILS_H_
#define UTILS_H_

#include "pyksp.h"
#include "globals.h"

void print_distances(std::vector<int> distances, const MyGraph & g);
Path pred_to_path(std::vector<std::size_t> preds, const MyGraph & g, Vertex source, Vertex sink);
void print_path(Path path, const MyGraph & g);

#endif
