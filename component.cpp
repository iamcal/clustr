/*
Clustr.  Copyright (c) 2007-2008 Yahoo! Inc.

All rights reserved.  This code is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License (GPL),
version 2 only.  This code is distributed WITHOUT ANY WARRANTY, whether
express or implied. See the GNU GPL for more details
(http://www.gnu.org/licenses/gpl.html)

AS A SPECIAL EXCEPTION, YOU HAVE PERMISSION TO LINK THIS PROGRAM WITH THE
CGAL LIBRARY AND DISTRIBUTE EXECUTABLES, AS LONG AS YOU FOLLOW THE REQUIREMENTS
OF THE GNU GPL VERSION 2 IN REGARD TO ALL OF THE SOFTWARE IN THE EXECUTABLE
ASIDE FROM CGAL. 
*/

#include <iostream>
#include "clustr.h"

using namespace Clustr;

Vertex_handle Mesh::vertex_circulator::target (const Edge &e) {
    return e.first->vertex(Mesh::ccw(e.second));
}

Vertex_handle Mesh::vertex_circulator::origin (const Edge &e) {
    return e.first->vertex(Mesh::cw(e.second));
}

Mesh::vertex_circulator&
Mesh::vertex_circulator::operator++(void)
{
    Mesh *mesh = ring->mesh;
    Edge_circulator edge = mesh->incident_edges(vertex), edge0 = edge;
    Edge selected;
    int min_count = ~0;

    if (previous.first == NULL) previous = *edge;
    Direction angle_in(mesh->segment(previous)), angle_out;

    (*vertex)++;
    do {
        if (mesh->classify(*edge) == Mesh::REGULAR) {
            Vertex_handle next_v = target(*edge);
            if (next_v == origin(previous)) continue;
            if (next_v->count() <= min_count) {
                Direction next_angle = mesh->segment(*edge).direction();
                if (next_v->count() == min_count && 
                    !next_angle.counterclockwise_in_between(angle_in,angle_out))
                    continue;
                selected  = *edge;
                angle_out = mesh->segment(selected).direction();
                min_count = next_v->count();
            }
        } 
    } while (++edge != edge0);
    //std::cerr << mesh->segment(previous) << " -> " << mesh->segment(selected) << std::endl;
    previous = selected;
    vertex = target(selected);
    return *this;
}
