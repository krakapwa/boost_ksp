#ifndef UTILS_H_
#define UTILS_H_

#include "pyksp.h"
#include "globals.h"

void print_distances(std::vector<int> distances, const MyGraph & g);
VertexPath pred_to_path(std::vector<std::size_t> preds, const MyGraph & g, Vertex source, Vertex sink);
void print_path(VertexPath path, const MyGraph & g);
void print_path(EdgePath path, const MyGraph & g);
void print_all(const MyGraph & g);
EdgePath vertpath_to_edgepath(VertexPath path, const MyGraph & g);

#endif
