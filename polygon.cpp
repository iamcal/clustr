/*
Clustr.  Copyright (c) 2007-2009 Yahoo! Inc.

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

# include <numeric>
# include <cmath>
# include "clustr.h"

using namespace Clustr;
//using std::tr1::sqrt;

void Ring::push_back (const Point &p2) {
    if (!is_empty()) {
        Point p1 = vertex(size()-1);
        coord_type fx = scale_x((p1.y() + p2.y())/2.0),
                   dx = p2.x() - p1.x(),
                   dy = p2.y() - p1.y();
        perimeter_ += sqrt(dx*dx*fx*fx + dy*dy);
        area_ += fx*fx*(p1.x() * p2.y() - p2.x() * p1.y()) / 2;
    }
    this->Ring_base::push_back(p2);
}

void Polygon::push_back (Ring &ring) {
    if ((empty() && !ring.is_ccw()) ||
        (!empty() && ring.is_ccw()))
        ring.reverse_orientation();
    this->Polygon_base::push_back(ring); 
}

coord_type Polygon::area (void) {
    coord_type total = 0;
    for (Polygon::iterator ring = begin(); ring != end(); ring++) {
        total += ring->area();
    }
    return total;
}

coord_type Polygon::perimeter (void) {
    coord_type total = 0;
    for (Polygon::iterator ring = begin(); ring != end(); ring++) {
        total += ring->perimeter();
    }
    return total;
}
