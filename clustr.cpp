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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#include <unistd.h> 
#include "clustr.h"
#include "shapefile.h"

using namespace std;
using namespace Clustr;

void write_points_to_shapefile (Shapefile &shape, vector<Point> &pts, 
                                string text) {
    vector<Point>::iterator begin = pts.begin(),
                            end   = pts.end();
    for (vector<Point>::iterator it = begin; it != end; it++) {
        Feature pt_feat(shape);
        Geometry pt_geom(wkbPoint);
        pt_geom.push_back(*it);
        pt_feat.set(pt_geom);
        pt_feat.set("tag", text.c_str());
        shape.add_feature(pt_feat);
    }
}

void write_polygon_to_shapefile (Shapefile &shape, Polygon &poly, string text, long int count) {
    Feature feat(shape);
    Geometry geom(wkbPolygon);
    geom.insert_rings(poly.begin(), poly.end());
    feat.set(geom);
    feat.set("tag", text.c_str());
    feat.set("count", count);
    feat.set("area", poly.area());
    feat.set("perimeter", poly.perimeter());
    feat.set("density", ((float)count)/poly.area());
    shape.add_feature(feat);
}

bool extract_alpha_component (Mesh::component &c, Polygon &poly, bool verbose)
{
    Mesh::vertex_circulator v(c);
    Vertex_handle v0(*v);
    Ring ring;

    do {
        ring.push_back(v->point());
    } while (*(++v) != v0);

    if (ring.size() < 4) {
        if (verbose) cerr << "Discarding degenerate shape." << endl;
        return false;
    } 
    if (verbose)
        cerr << "- " << ring.size() << " vertices, "
             << "area: " << fabs(ring.area())
             << ", perimeter: " << ring.perimeter() << endl;

    poly.push_back(ring);
    return true;
}

bool compare_polygon_sizes (const Polygon &a, const Polygon &b) {
    return a.size() < b.size();
}

void attach_interior_rings (vector<Polygon> &input, vector<Polygon> &output, bool verbose = false) {
    unsigned int count = 0;

    /* First, sort the polygons by area. This guarantees that smaller polygons
     * will wind up inside larger polygons. */
    sort(input.begin(), input.end(), compare_polygon_sizes);

    for (vector<Polygon>::iterator it = input.begin(); it != input.end(); it++) {
        bool is_interior_ring = false;

        for (vector<Polygon>::iterator comp = it + 1; comp != input.end(); comp++) {
            /* check if the first point in the first ring of *it is inside the outer ring
             * of *comp. If it is, then we assume that the entire polygon *it is actually
             * inside *comp.
             *
             * Note that *comp needs to be a simple polygon -- it *should* be but sometimes
             * the polygons that come out of the alpha shape self-intersect at a point.
             * Not sure what to do about this. */
            if ((*comp)[0].is_simple() && (*comp)[0].encloses((*it)[0][0])) {
                for (Polygon::iterator ring = it->begin(); ring != it->end(); ring++) {
                    // this is actually going to get double-interior rings backwards,
                    // but fixing that means making the logic in
                    // Polygon::push_back a lot more sophisticated...
                    comp->push_back(*ring); 
                }
                is_interior_ring = true;
                count++;
                break;
            }
        }

        if (!is_interior_ring)
            output.push_back(*it);
    }

    if (verbose && count)
        cerr << count << " interior rings attached." << endl;
}

void construct_alpha_shape (Config &config, Shapefile &shape,
                            vector<Point> &pts, string text)
{
    coord_type alpha = config.alpha;
    Mesh mesh(pts.begin(), pts.end(), alpha);
    vector<Polygon> polygons, output;

    if (alpha == 0) {
        if (config.verbose) cerr << "Computing optimal alpha... ";
        alpha = *mesh.find_optimal_alpha(1);
        if (config.verbose) cerr << alpha << endl;
        mesh.set_alpha(alpha);
    }

    if (config.verbose)
        cerr << mesh.number_of_solid_components() <<
              " component(s) found for alpha value " << alpha << "." << endl;

    for (Mesh::component c = mesh.components_begin(); c != mesh.components_end(); c++) {
        Polygon poly;
        if (extract_alpha_component(c, poly, config.verbose))
            polygons.push_back(poly);
    }

    attach_interior_rings(polygons, output, config.verbose);

    if (config.verbose) cerr << "Writing " << output.size() << " polygons to shapefile." << endl;
    for (vector<Polygon>::iterator poly = output.begin(); poly != output.end(); poly++) {
        write_polygon_to_shapefile(shape, *poly, text, pts.size()); // FIXME: size is broken
    }
}

void construct_output (Config &config, Shapefile &shape,
                       vector<Point> &pts, string tag) {
    cerr << "Got " << pts.size() << " points for tag '"
                   << tag << "'." << endl;
    if (config.points_only) {
        write_points_to_shapefile(shape, pts, tag);
    } else {
        if (pts.size() >= 3) {
            construct_alpha_shape(config, shape, pts, tag);
        } else {
            if (config.verbose)
                cerr << "Need more than 3 points for tag '" 
                     << tag << "'" << endl;
        }
    }
}

void display_usage (void) {
    cerr << endl
         << "clustr "  << VERSION << " - construct polygons from tagged points"
         << endl
         << "written by " << AUTHOR << endl
         << COPYRIGHT  << endl
         << endl
         << "Usage: clustr [-a <n>] [-p] [-v] <input> <output>" << endl
         << "   -h, -?      this help message" << endl
         << "   -v          be verbose (default: off)" << endl
         << "   -a <n>      set alpha value (default: use \"optimal\" value)"
         << endl
         << "   -p          output points to shapefile, instead of polygons" 
         << endl 
         << endl
         << "If <input> is missing or given as \"-\", stdin is used."
         << endl
         << "If <output> is missing, output is written to clustr.shp."
         << endl
         << "Input file should be formatted as: <tag> <lon> <lat>\\n" << endl
         << "Tags must not contain spaces." << endl
         << endl;
}

bool parse_options (Config &config, int argc, char **argv) {
    istringstream iss;
    string optstring("a:pvh?");

    if (argc <= 1) {
        display_usage();
        return false;
    }

    int opt = getopt( argc, argv, optstring.c_str() );
    while( opt != -1 ) {
        switch( opt ) {
            case 'a':
                iss.str((string) optarg);
                iss >> config.alpha;
                break;
            case 'p':
                config.points_only = true;
                break;
            case 'v':
                config.verbose = true;
                break;
            case 'h':   /* fall-through is intentional */
            case '?':
                display_usage();
                return false;
        }
        opt = getopt( argc, argv, optstring.c_str() );
    }
    if (optind < argc) {
        config.in_file = argv[optind++];
    } 
    if (optind < argc) {
        config.out_file = argv[optind];
    } 
    return true;
}

int main(int argc, char **argv) {
    vector<Point> points;
    istream *input;
    Config config;
    Shapefile *shape;

    if (!parse_options(config, argc, argv))
        return 1;

    if (config.in_file == "-") {
        input = &cin;
    } else {
        input = new ifstream(config.in_file.c_str());
    }
    
    if (config.points_only) {
        if (config.verbose) {
            cerr << "Only converting points to Shapefile,"
                 << " no clustering performed." << endl;
        }
        shape = new Shapefile(config.out_file, wkbPoint);
        shape->add_field("tag", OFTString, 64);
        
    } else {
        shape = new Shapefile(config.out_file, wkbPolygon);
        shape->add_field("tag", OFTString, 64);
        shape->add_field("count", OFTInteger, 8);
        shape->add_field("area", OFTReal, 10, 2);
        shape->add_field("perimeter", OFTReal, 10, 2);
        shape->add_field("density", OFTReal, 10, 2);
    }

    string line, tag, previous = "";
    vector<Point> pts;
    Point pt;
    
    if (config.verbose)
        cerr << "Reading points from input." << endl;

    while (!getline(*input,line).eof()) {
        istringstream is(line);
        if (!(is >> tag >> pt)) {
            cerr << "Bad input line: " << line;
            continue;
        }
        if (previous == "") {
            previous = tag;
        }
        if (tag != previous) {
            construct_output(config, *shape, pts, previous);
            previous = tag;
            pts.clear();
        }
        pts.push_back(pt);
    }
    construct_output(config, *shape, pts, tag);
   
    // if you want to do memory leak testing, include <ogr_api.h> and 
    // uncomment the following line:
    // OGRCleanupAll();
    delete shape;
    return 0;
}

