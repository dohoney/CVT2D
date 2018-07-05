/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#ifndef __DELAUNAY_CGAL__
#define __DELAUNAY_CGAL__

#include "geometry.h"

//______________________________________________________________________________________
// CGAL stuff

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
///dxy add
#include <CGAL/Polygon_2.h>
#include <CGAL/Object.h>
///

namespace Geex {


// ------------- My stuff to extend CGAL --------------------------

    class VertexDecoration {
    public:
        int index ;
        bool dual_intersects_boundary ;
        double theta ;
        double rho ;
        double energy ;
        bool locked ;
		///dxy add
		vec2 lloyd_grad ;
		vec2 direction_grad;
		double area;
		bool selected;
		//int valence;
		int singularity; //used to mark different types of singularity
		int separatrix_in_degree; //used to compute intersections of separatrix
		///
    } ;

    class CellDecoration {
    public:

        CellDecoration() {
        }

        bool infinite ;
        vec2 dual ;
        bool dual_outside ;
		///dxy add
		bool selected;
		bool selected_edge[3];
		//edge stretch
		double edge_stretch[3];
		int group_mark;
		///
    } ;


// ------------- CGAL stuff ---------------------------------

    // ----------------------- A CGAL::Vertex with decoration ------------------
    template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
    class Vertex : public  Vb, public VertexDecoration  {
        typedef Vb superclass;
    public:
        typedef typename Vb::Vertex_handle      Vertex_handle;
        typedef typename Vb::Face_handle        Face_handle;
        typedef typename Vb::Point              Point;
        
        template < typename TDS2 >
        struct Rebind_TDS {
            typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
            typedef Vertex<Gt,Vb2> Other;
        } ;
        
    public:
        Vertex() : superclass() {}
        Vertex(const Point & p) : superclass(p) {}
        Vertex(const Point & p, Face_handle f) : superclass(f,p) {}
        Vertex(Face_handle f) : superclass(f) {}
    } ;


    // ----------------------- A CGAL::Cell with decoration ------------------
    template < class Gt, class Cb = CGAL::Triangulation_face_base_2<Gt> >
    class Cell : public Cb, public CellDecoration {
        typedef Cb superclass;
    public:
        typedef typename Cb::Vertex_handle      Vertex_handle;
        typedef typename Cb::Face_handle        Face_handle;
        template < typename TDS2 >
        struct Rebind_TDS {
            typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
            typedef Cell<Gt,Cb2> Other;
        };


        Cell() : superclass() {  }
        Cell(
            Vertex_handle v0, Vertex_handle v1, Vertex_handle v2
        ) : superclass(v0,v1,v2) { }
            
        Cell(
            Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, 
            Face_handle n0, Face_handle n1, Face_handle n2 
        ) : superclass(v0,v1,v2,n0,n1,n2) { }

    } ;
    
    typedef double Coord_type;
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;
    typedef Vertex<K> Vb ;
    typedef  CGAL::Triangulation_hierarchy_vertex_base_2<Vb> Vbh;
    typedef Cell<K>   Cb ;
    typedef CGAL::Triangulation_data_structure_2<Vbh,Cb> TDS;
    typedef CGAL::Delaunay_triangulation_2<K, TDS> TRI ;

    typedef K::Point_2    Point;
    typedef K::Vector_2   CGAL_Vector;
	///dxy add
	typedef K::Segment_2	Segment_2;
	typedef K::Ray_2		Ray_2;
	typedef CGAL::Polygon_2<K> Polygon_2;
	///

    inline vec2 to_geex(const Point& p) { return vec2(p.x(), p.y()) ; }
    inline Point to_cgal(const vec2& p) { return Point(p.x, p.y) ; }

    class DelaunayBase : public CGAL::Triangulation_hierarchy_2<TRI> {
    public:
    } ;

}


#endif
