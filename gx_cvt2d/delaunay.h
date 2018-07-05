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

#ifndef __DELAUNAY__
#define __DELAUNAY__

//#define NOMINMAX 

#include "delaunay_cgal.h"
#include "polygons.h"
#include <GL/gl.h>
///dxy add
#include <queue>
///

namespace Geex {

	class DelaunayGraphics ;
	class DelaunayCVT ;
	//______________________________________________________________________________________

	class SegmentDelaunay : public DelaunayBase {
		typedef DelaunayBase baseclass ;
	public:
		void insert(int id, const vec2& p1, const vec2& p2, double step_length) ;
		int locate(const vec2& p) ;
	protected:
		void insert(int id, const vec2& p) ;
	} ;

	class Delaunay : public DelaunayBase {
		typedef DelaunayBase baseclass ;

	public:
		Delaunay() ;

		void save(const std::string& filename) ;
		void load(const std::string& filename) ;
		void load_boundary(const std::string& filename) ;
		//dxy add
		void load_boundary(const std::string& filename, int nb_pts);
		//
		void set_non_convex_mode(bool x) { non_convex_mode_ = x ; }

		void get_bbox(
			real& x_min, real& y_min, real& z_min,
			real& x_max, real& y_max, real& z_max
			) ;

		bool in_boundary(const vec2& p) { 
			return non_convex_mode_ ? boundary_.contains(p) : boundary_convex_.contains(p) ; 
		}

		GLboolean& insert_boundary() { return insert_boundary_ ; }

		// ------------------------------------ Delaunay 

		int nb_vertices() const { return all_vertices_.size() ; }
		void clear() ;
		void begin_insert() ;
		Vertex_handle insert(const vec2& p) ;
		Vertex_handle nearest(const vec2& p) ;
		void remove(const vec2& p) ;
		void end_insert(bool redraw = true) ;

		void insert_random_vertex() ;
		void insert_random_vertices(int nb) ;
		void insert_grid() ;

		std::vector<Vertex_handle>& vertices() { return all_vertices_ ; }

		vec2 vertex(unsigned int i) const {
			return vec2(
				all_vertices_[i]->point().x(),
				all_vertices_[i]->point().y()
				) ;
		}

		// ------------------------------------ Delaunay combinatorics ----------------

		bool dual_cell_intersects_boundary(Vertex_handle v) {
			return v->dual_intersects_boundary ;
		}

		Line<real> get_dual_line(Edge e) {
			return median_line(
				to_geex(e.first->vertex(ccw(e.second))->point()),
				to_geex(e.first->vertex( cw(e.second))->point())
				) ;
		}
		Polygon2* dual_convex_clip(Vertex_handle v, bool close = true) ;
		int dual_facet_degree(Vertex_handle v) ;

		// --------------

		SegmentDelaunay& segments() { return segments_ ; }

		///dxy add
		Polygon_2 to_cgal_polygon(const Polygon2& png);
		double boundary_clipped_length(const CGAL::Object& obj);
		
		//short long topo edit
		void split_long_edge();
		void collapse_short_edge();
		void split_isolate_long_edge();
		void collapse_isolate_short_edge();
		//select
		enum SelectMode {VERTEX, EDGE, FACE, SEPARATRIX_SEGMENT };
		GLenum& select_mode() { return select_mode_ ; }
		void clear_select();
		void update_select(const vec2&);
		//selected edge edit
		void split_selected_edge();
		void collapse_selected_edge();
		void flip_selected_edge();
		//stretch topo edit
		void calc_edge_stretch();
		void calc_vertex_singularity();
		void calc_face_group_mark();  // faces into non-adjacent groups
		void stretch_topo_optimize();
		/*bool& is_edge_stretch_dirty() { return is_edge_stretch_dirty_; }
		bool& is_vertex_singularity_dirty() { return is_vertex_singularity_dirty_; }
		bool& is_face_group_mark_dirty() { return is_face_group_mark_dirty_; }*/
		int face_group_size() { 
			calc_face_group_mark();
			return face_group_size_; 
		}
		bool has_irregular(Delaunay::Face_handle);
		//
		void update_area(bool info = true);
		///dxy add end

	protected:
		std::vector<Vertex_handle> all_vertices_ ;
		bool non_convex_mode_ ;

		Polygon2 boundary_ ;
		Convex boundary_convex_ ;
		bool cached_bbox_ ;
		double x_min_, y_min_, x_max_, y_max_ ;

		Polygon2 ping_ ;
		Polygon2 pong_ ;

		SegmentDelaunay segments_ ;

		friend class DelaunayGraphics ;
		friend class DelaunayCVT ;

		bool opened_ ;
		GLboolean insert_boundary_ ;

		///dxy add
		Polygon_2 cgal_boundary_ ;

		bool is_edge_stretch_dirty_;
		bool is_vertex_singularity_dirty_;
		bool is_face_group_mark_dirty_;
		int face_group_size_;

		float total_area_; //used to normalize cvt and direction energy
		float average_area_;

		GLenum select_mode_;
		///
	} ;

	//______________________________________________________________________________________

#define FOR_EACH_VERTEX(TCLASS, INSTANCE, IT)                       \
	for(                                                            \
	TCLASS::Vertex_iterator IT = (INSTANCE)->vertices_begin() ; \
	IT != (INSTANCE)->vertices_end(); IT++                      \
	)

#define FOR_EACH_FACE(TCLASS, INSTANCE, IT)                      \
	for(                                                         \
	TCLASS::Face_iterator IT = (INSTANCE)->faces_begin() ;   \
	IT != (INSTANCE)->faces_end(); IT++                      \
	)


	///dxy add: triEdge 
	class triEdge {
	public:
		triEdge(Delaunay::Face_handle f = nullptr, unsigned idx = 0) : face_(f), edge_index_(idx) {}
		Delaunay::Vertex_handle start() const { return face_->vertex((edge_index_ + 1) % 3); }
		Delaunay::Vertex_handle end() const { return face_->vertex((edge_index_ + 2) % 3); }
		Delaunay::Face_handle face() const { return face_; }
		unsigned edge_index() const { return edge_index_; }
		bool operator==(const triEdge &rhs) const {
			return (this->start() == rhs.start() && this->end() == rhs.end())
				|| (this->start() == rhs.end() && this->end() == rhs.start());
		}
		bool operator!=(const triEdge &rhs) const {	return !(*this == rhs);	}
		bool is_valid() const { return (face_ != nullptr) && (edge_index_ <= 3); }
	private:
		Delaunay::Face_handle face_;
		unsigned edge_index_;
	};

	typedef std::pair<triEdge, double> EdgeValue;
	struct EdgeValue_greater{  
		bool operator ()(EdgeValue &a, EdgeValue &b){  
			return a.second > b.second;  
		}  
	}; 
	struct EdgeValue_less{  
		bool operator ()(EdgeValue &a, EdgeValue &b){  
			return a.second < b.second;  
		}  
	};
	///
}

#endif
