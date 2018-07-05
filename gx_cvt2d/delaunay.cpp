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
*  You should have received a copy of the GNU General Public License<
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

#include "delaunay.h"
#include <Geex/basics/file_system.h>
#include <glut_viewer/glut_viewer.h>
#include <fstream>

namespace Geex {

	void SegmentDelaunay::insert(int id, const vec2& p1, const vec2& p2, double step_length) {
		vec2 V = p2 - p1 ;
		double l = V.length() ;
		int count = int(l / step_length) ;
		count = gx_max(count, 2) ;
		V = normalize(V) ;
		double delta = l / double(count) ;
		for(int i=0; i<=count; i++) {
			insert(id, p1 + (delta * double(i)) * V) ;
		}
	} 

	int SegmentDelaunay::locate(const vec2& p) {
		return nearest_vertex(to_cgal(p))->index ;
	}

	void SegmentDelaunay::insert(int id, const vec2& p) {
		baseclass::insert(to_cgal(p))->index = id ;
	}

	//____________________________________________________________________________________

	Delaunay::Delaunay() {
		non_convex_mode_ = false ;
		cached_bbox_ = false ;
		load_boundary(FileSystem::get_project_root() + "/gx_cvt2d/square.line") ;
		opened_ = false ;
		insert_boundary_ = false ;
	}

	void Delaunay::save(const std::string& filename) {
		std::ofstream out(filename.c_str()) ;
		for(Vertex_iterator it = vertices_begin() ; it != vertices_end() ; it++) {
			if(is_infinite(it)) { continue ; }
			out << to_geex(it->point()) << std::endl ;
		}
	}

	void Delaunay::load(const std::string& filename) {
		std::cerr << "loading " << filename << std::endl ;
		clear() ;
		std::ifstream in(filename.c_str()) ;
		if(!in) {
			std::cerr << "could not open file" << std::endl ;
			return ;
		}
		vec2 p ;
		while(in) {
			in >> p ;
			if(in) {  // we need to do this additional check else we got the last point twice !
				insert(p) ; 
			}
		} ;
	}

	void Delaunay::load_boundary(const std::string& filename) {
		cached_bbox_ = false ;
		boundary_.clear() ;
		boundary_convex_.clear() ;
		boundary_.load(filename) ;
		//boundary_.normalize() ;
		if(!non_convex_mode_) {
			for(unsigned int i=0; i<boundary_.size(); i++) {
				boundary_convex_.push_back(boundary_[i].line()) ;
			}
		}

		double x_min, y_min, z_min ;
		double x_max, y_max, z_max ;
		get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
		double dx = x_max - x_min ;
		double dy = y_max - y_min ;
		double diag = sqrt(dx*dx + dy*dy) ;

		segments_.clear() ;
		for(unsigned int i=0; i<boundary_.size(); i++) {
			segments_.insert(int(i), boundary_[i].vertex[0], boundary_[i].vertex[1], 0.01 * diag) ;
		}

		///dxy add
		cgal_boundary_.clear();
		cgal_boundary_ = to_cgal_polygon(boundary_);
		///

	}

	///dxy add
	void Delaunay::load_boundary(const std::string& filename, int nb_pts) {
		cached_bbox_ = false ;
		boundary_.clear() ;
		boundary_convex_.clear() ;
		boundary_.load(filename) ;
		//dxy add
		boundary_.cell_normalize(nb_pts) ;
		//
		if(!non_convex_mode_) {
			for(unsigned int i=0; i<boundary_.size(); i++) {
				boundary_convex_.push_back(boundary_[i].line()) ;
			}
		}

		double x_min, y_min, z_min ;
		double x_max, y_max, z_max ;
		get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
		double dx = x_max - x_min ;
		double dy = y_max - y_min ;
		double diag = sqrt(dx*dx + dy*dy) ;

		segments_.clear() ;
		for(unsigned int i=0; i<boundary_.size(); i++) {
			segments_.insert(int(i), boundary_[i].vertex[0], boundary_[i].vertex[1], 0.01 * diag) ;
		}

		///dxy add
		cgal_boundary_.clear();
		cgal_boundary_ = to_cgal_polygon(boundary_);
		///
	}
	///dxy add end

	void Delaunay::get_bbox(
		real& x_min, real& y_min, real& z_min,
		real& x_max, real& y_max, real& z_max
		) {
			z_min = 0.0 ; z_max = 1.0 ;
			if(!cached_bbox_) {
				x_min_ = y_min_ =  1e30 ;
				x_max_ = y_max_ = -1e30 ;
				if(boundary_.size() == 0) {
					x_min_ = 0.0 ; y_min_ = 0.0 ;
					x_max_ = 1.0 ; y_max_ = 1.0 ;
				}
				for(unsigned int i=0; i<boundary_.size(); i++) {
					for(unsigned int j=0; j<2; j++) {
						x_min_ = gx_min(x_min_, boundary_[i].vertex[j].x) ;
						y_min_ = gx_min(y_min_, boundary_[i].vertex[j].y) ;
						x_max_ = gx_max(x_max_, boundary_[i].vertex[j].x) ;
						y_max_ = gx_max(y_max_, boundary_[i].vertex[j].y) ;
					}
				}
				cached_bbox_ = true ;
			}
			x_min = x_min_ ; y_min = y_min_ ;
			x_max = x_max_ ; y_max = y_max_ ;
	}

	// ------------------------------------ Delaunay 

	void Delaunay::clear() {
		baseclass::clear() ;
		all_vertices_.clear() ;
	}

	Delaunay::Vertex_handle Delaunay::insert(const vec2& p) {
		gx_assert(opened_) ;
		if(Numeric::is_nan(p.x) || Numeric::is_nan(p.y)) {
			std::cerr << "Nan !" << std::endl ;
			return 0 ;
		}
		Vertex_handle result = baseclass::insert(to_cgal(p)) ;
		result->theta = ::Geex::Numeric::random_float64() * M_PI / 2.0 ;
		result->energy = 0.0 ;
		result->rho = 1.0 ; 
		result->locked = false ;
		///dxy add
		result->area = 0.0;
		///
		return result ;
	}

	Delaunay::Vertex_handle Delaunay::nearest(const vec2& p) {
		Delaunay::Face_handle f = locate(to_cgal(p)) ;
		Delaunay::Vertex_handle result = f->vertex(0) ;
		double dist = (to_geex(result->point()) - p).length2() ;
		for(unsigned int i=1; i<3; i++) {
			double cur_dist = (to_geex(f->vertex(i)->point()) - p).length2() ;
			if(cur_dist < dist) {
				dist = cur_dist ;
				result = f->vertex(i) ;
			}
		}
		return result ;
	}

	void Delaunay::remove(const vec2& p) {
		if(number_of_vertices() <= 3) { return ; }
		gx_assert(opened_) ;
		Face_handle f = locate(to_cgal(p)) ;


		double min_d = 1e30 ;
		Vertex_handle v = 0 ;
		for(unsigned int i=0; i<3; i++) {
			if(!is_infinite(f->vertex(i))) {
				double cur_d = (to_geex(f->vertex(i)->point()) - p).length2() ;
				if(cur_d < min_d) {
					min_d = cur_d ;
					v = f->vertex(i) ;
				}
			}
		}
		baseclass::remove(v) ;

	}

	void Delaunay::begin_insert() { opened_ = true ; }

	///dxy change
	void Delaunay::end_insert(bool redraw) {
		all_vertices_.clear() ;
		for(All_vertices_iterator it = all_vertices_begin() ; it != all_vertices_end() ; it++) {
			it->dual_intersects_boundary = false ;
			it->index = -1 ;
			if(!is_infinite(it)) {
				all_vertices_.push_back(it) ;
			}
		}
		for(All_faces_iterator it = all_faces_begin() ; it != all_faces_end() ; it++) {
			if(is_infinite(it)) {
				it->infinite = true ;
				it->dual = vec2(0.0, 0.0) ;
				it->dual_outside = true ;
			} else {
				it->infinite = false ;
				it->dual = to_geex(baseclass::dual(it)) ;
				it->dual_outside = !in_boundary(it->dual) ;
			}
			if(it->dual_outside) {
				it->vertex(0)->dual_intersects_boundary = true ;
				it->vertex(1)->dual_intersects_boundary = true ;
				it->vertex(2)->dual_intersects_boundary = true ;
			}
		}        

		all_vertices_.clear() ;
		int cur_index = 0 ;
		for(Vertex_iterator it = vertices_begin(); it != vertices_end() ; it++) {
			all_vertices_.push_back(it) ;
			it->index = cur_index ;
			cur_index++ ;
		}
		opened_ = false ;
		if(redraw) {
			glut_viewer_redraw() ;            
		}

		///dxy add: set states
		clear_select();
		is_edge_stretch_dirty_ = true;
		is_vertex_singularity_dirty_ = true;
		is_face_group_mark_dirty_ = true;
		face_group_size_ = 0;
		///
	}


	inline vec2 random_v() {
		return vec2( 
			Numeric::random_float64(),
			Numeric::random_float64()
			) ;
	}


	void Delaunay::insert_random_vertex() {
		double x_min, y_min, z_min, x_max, y_max, z_max ;
		get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
		Geex::vec2 p = random_v() ;
		p.x = x_min + (x_max - x_min) * p.x ;
		p.y = y_min + (y_max - y_min) * p.y ;
		int nb_tries = 0 ;
		while(!in_boundary(p)) {
			nb_tries++ ;
			if(!non_convex_mode_ && nb_tries > 1000) {
				std::cerr << "Could not insert point, probably missing +non_convex flag" << std::endl ;
				exit(-1) ;
			}
			p = random_v() ;
			p.x = x_min + (x_max - x_min) * p.x ;
			p.y = y_min + (y_max - y_min) * p.y ;
		}
		insert(p) ;
	}

	void Delaunay::insert_random_vertices(int nb) {
		begin_insert() ;
		for(unsigned int i=0; i<nb; i++) {
			insert_random_vertex() ;
			//            std::cerr << (i+1) << '/' << nb << std::endl ;
		}
		if(insert_boundary_ ) {
			for(unsigned int i=0; i<boundary_.size(); i++) {
				Vertex_handle v = insert(boundary_[i].vertex[0]) ;
				v->locked = true ;
			}
		}
		end_insert(false) ;
	}

	void Delaunay::insert_grid() {

		begin_insert() ;
		for(unsigned int i=0; i<10; i++) {
			for(unsigned int j=0; j<10; j++) {
				double x = 0.05 + double(i) / double(10.0) ;
				double y = 0.05 + double(j) / double(10.0) ;
				insert(vec2(x,y)) ;
			}
		}
		end_insert(false) ;

		/*
		double x_min, y_min, z_min, x_max, y_max, z_max ;
		get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
		double d = x_max - x_min ;
		x_min -= 0.1*d ;
		x_max += 0.1*d ;
		y_min -= 0.1*d ;
		y_max += 0.1*d ;
		begin_insert() ;
		for(unsigned int i=0; i<16; i++) {
		for(unsigned int j=0; j<16; j++) {
		double x = x_min + double(i) / double(15.0) * (x_max - x_min) ;
		double y = y_min + double(j) / double(15.0) * (y_max - y_min) ;
		insert(vec2(x,y)) ;
		}
		}
		end_insert(false) ;
		*/
	}


	Polygon2* Delaunay::dual_convex_clip(Vertex_handle v, bool close) {
		Polygon2* from = &boundary_ ;
		Polygon2* to = &pong_ ;

		std::vector<Edge> edges ;
		Edge_circulator it = incident_edges(v) ;
		do {
			edges.push_back(*it) ;
			it++ ;
		} while(it != incident_edges(v)) ;

		// No need to have special cases for infinite vertices: they
		// do not yield any clipping plane !
		for(unsigned int i=0; i<edges.size(); i++) {
			Edge e = edges[i] ;
			if(is_infinite(e)) { continue ; }
			Geex::Line<real> L = get_dual_line(e) ;
			int E = e.first->vertex(ccw(e.second))->index ;
			gx_assert(E != v->index) ;
			from->convex_clip(*to, L, E, close) ;
			if(from == &boundary_) {
				from = &pong_ ;
				to = &ping_ ;
			} else {
				gx_swap(from, to) ;
			}
		}
		return from ;
	}

	int Delaunay::dual_facet_degree(Vertex_handle v) {
		int result = 0 ;
		if(dual_cell_intersects_boundary(v)) { 
			Polygon2* P = dual_convex_clip(v) ;
			result = P->size() ;
		} else {
			Face_circulator it = incident_faces(v) ;
			do {
				result++ ;
				it++ ;
			} while(it != incident_faces(v)) ;
		} 
		return result ;
	}

	//////////////////////////////////////////////////////////////////////////
	///dxy add
	void Delaunay::split_long_edge() {
		std::set<Delaunay::Vertex_handle> edgeStart;
		std::vector<vec2> to_insert;
		FOR_EACH_FACE(Delaunay, this, it) {
			if (is_infinite(it)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = i;
				int j2 = (i+1) % 3;
				Delaunay::Vertex_handle v1 = it->vertex(j1);
				Delaunay::Vertex_handle v2 = it->vertex(j2);
				if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue;  //exclude boundary vertex
				int d1 = dual_facet_degree(v1);
				int d2 = dual_facet_degree(v2);

				//if (d1 > 6 && d2 > 6)  //long edge condition 0
				if (d1 + d2 >= 14) //long edge condition 1
				{
					if (edgeStart.find(v2) == edgeStart.end())
					{
						edgeStart.insert(v1);
						double xnew = 0.5 * (v1->point().x() + v2->point().x());
						double ynew = 0.5 * (v1->point().y() + v2->point().y());
						to_insert.push_back(vec2(xnew, ynew));
					}
				}
			}
		}
		//insert
		begin_insert() ;
		for (int i=0; i<to_insert.size(); ++i)
		{
			insert(to_insert[i]);
		}
		end_insert(false) ;
		std::cout << to_insert.size() << " points inserted." << std::endl;
	}

	void Delaunay::split_isolate_long_edge() {
		std::set<Delaunay::Vertex_handle> incident_vertices;
		std::vector<vec2> to_insert;
		FOR_EACH_FACE(Delaunay, this, it) {
			if (is_infinite(it)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = i;
				int j2 = (i+1) % 3;
				Delaunay::Vertex_handle v1 = it->vertex(j1);
				Delaunay::Vertex_handle v2 = it->vertex(j2);
				if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue;  //exclude boundary vertex
				int d1 = dual_facet_degree(v1);
				int d2 = dual_facet_degree(v2);

				//if (d1 > 6 && d2 > 6)  //long edge condition 0
				if (d1 + d2 >= 14) //long edge condition 1
				{
					if (incident_vertices.find(v1) == incident_vertices.end() &&
						incident_vertices.find(v2) == incident_vertices.end()) //isolate condition
					{
						incident_vertices.insert(v1);
						incident_vertices.insert(v2);
						double xnew = 0.5 * (v1->point().x() + v2->point().x());
						double ynew = 0.5 * (v1->point().y() + v2->point().y());
						to_insert.push_back(vec2(xnew, ynew));
					}
				}
			}
		}
		//insert
		begin_insert() ;
		for (int i=0; i<to_insert.size(); ++i)
		{
			insert(to_insert[i]);
		}
		end_insert(false) ;
		std::cout << to_insert.size() << " points inserted." << std::endl;
	}

	void Delaunay::collapse_short_edge() {
		std::set<Delaunay::Vertex_handle> to_remove;
		std::vector<vec2> to_insert;
		FOR_EACH_FACE(Delaunay, this, it) {
			if (is_infinite(it)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = i;
				int j2 = (i+1) % 3;
				Delaunay::Vertex_handle v1 = it->vertex(j1);
				Delaunay::Vertex_handle v2 = it->vertex(j2);
				if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue;  //exclude boundary vertex
				int d1 = dual_facet_degree(v1);
				int d2 = dual_facet_degree(v2);

				//if (d1 < 6 && d2 < 6)  //short edge condition 0
				if (d1 + d2 <= 10) //short edge condition 1
				{
					if (to_remove.find(v1) == to_remove.end() || to_remove.find(v2) == to_remove.end())
					{
						to_remove.insert(v1);
						to_remove.insert(v2);
						double xnew = 0.5 * (v1->point().x() + v2->point().x());
						double ynew = 0.5 * (v1->point().y() + v2->point().y());
						to_insert.push_back(vec2(xnew, ynew));
					}
				}
			}
		}
		//insert and remove
		begin_insert() ;
		for (int i=0; i<to_insert.size(); ++i)
		{
			insert(to_insert[i]);
		}
		for (auto it=to_remove.begin(); it != to_remove.end(); ++it)
		{
			remove(vec2((*it)->point().x(), (*it)->point().y()));
		}
		end_insert(false) ;
		std::cout << to_remove.size() - to_insert.size() << " points removed." << std::endl;
	}

	void Delaunay::collapse_isolate_short_edge() {
		collapse_short_edge(); //this already satisfies isolate condition
	}

	//predicates
	bool Delaunay::has_irregular(Delaunay::Face_handle f) { //user should guarantee vertex->singularity is up to date
		for (int i=0; i<3; ++i)
		{
			Delaunay::Vertex_handle v = f->vertex(i);

			//if (this->dual_facet_degree(v) != 6  && !v->dual_intersects_boundary)
			if (v->singularity != 0)
			{
				return true;
			}
		}
		return false;
	}

	//calc vertex singularity type
	void Delaunay::calc_vertex_singularity() {
		if (!is_vertex_singularity_dirty_)
		{
			return;
		}

		FOR_EACH_VERTEX(Delaunay, this, it) {
			if (this->dual_facet_degree(it) != 6  && !it->dual_intersects_boundary)  //inner singularity
			{
				it->singularity = 1;
			} 
			else
			{
				it->singularity = 0;
			}
		}

		/*FOR_EACH_FACE(Delaunay, this, it) {
			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v1 = it->vertex(i);
				Delaunay::Vertex_handle v2 = it->vertex((i+1)%3);
				if (v2->singularity == 0
					&& v1->singularity == 1
					&& v2->dual_intersects_boundary
					&& (this->dual_facet_degree(v2)+2) != 6)
				{
					v2->singularity = 2;
				}
			}
		}*/

		is_vertex_singularity_dirty_ = false;
	}

	//calc edge stretch
	void Delaunay::calc_edge_stretch() {
		if (!is_edge_stretch_dirty_)
		{
			return;
		}

		FOR_EACH_FACE(Delaunay, this, it) {
			if (this->is_infinite(it)) {
				it->edge_stretch[0] = 0;
				it->edge_stretch[1] = 0;
				it->edge_stretch[2] = 0;
			}

			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v1 = it->vertex((i+1)%3);
				Delaunay::Vertex_handle v2 = it->vertex((i+2)%3);
				vec2 p1 = to_geex(v1->point());
				vec2 p2 = to_geex(v2->point());
				vec2 v12 = p2 - p1;
				vec2 e12 = normalize(v12);

				it->edge_stretch[i] = dot(v1->direction_grad, e12) + dot(v2->direction_grad, -e12);
			}
		}
		is_edge_stretch_dirty_ = false;
	}

	//calc face group mark
	void Delaunay::calc_face_group_mark() {
		if (!is_face_group_mark_dirty_)
		{
			return;
		}
		calc_vertex_singularity(); //make sure vertex singularity is up to date

		//clear face group mark
		FOR_EACH_FACE(Delaunay, this, it) {
			it->group_mark = -1;
		}

		//calc face group mark
		int hmark = 0;
		FOR_EACH_FACE(Delaunay, this, it) {
			if (is_infinite(it)) continue;

			if (it->group_mark == -1 && has_irregular(it))
			{
				it->group_mark = hmark;
				//faces bfs
				std::queue<Delaunay::Face_handle> Q;
				Q.push(it);
				while (!Q.empty())
				{
					Delaunay::Face_handle t = Q.front();
					Q.pop();
					//visit and mark neighbor faces of t
					for (int i=0; i<3; ++i)
					{
						Delaunay::Vertex_handle v = t->vertex(i);
						Delaunay::Face_circulator s = incident_faces(v, t);
						++s;
						do 
						{
							if (s->group_mark == -1 && has_irregular(s))
							{
								s->group_mark = hmark;
								Q.push(s);
							}
							++s;
						} while (s != t);
					}
				}
				++hmark;
			}
		}
		face_group_size_ = hmark;
		is_face_group_mark_dirty_ = false;
		//std::cout << "hmark = " << hmark << std::endl;
	}

	//stretch topo optimize
	void Delaunay::stretch_topo_optimize() {
		calc_edge_stretch();
		calc_vertex_singularity();
		calc_face_group_mark();
		int hmark = face_group_size_;

		//calc stretched edge for each face group
		std::vector<EdgeValue> isolate_edges(hmark, std::make_pair(triEdge(), 0));

		FOR_EACH_FACE(Delaunay, this, it) {
			int gmark = it->group_mark;
			if (gmark != -1)
			{
				for (int i=0; i<3; ++i)
				{
					Delaunay::Vertex_handle v1 = it->vertex((i+1)%3);
					Delaunay::Vertex_handle v2 = it->vertex((i+2)%3);

					//
					if (v1->singularity != 0 
						|| (v1->dual_intersects_boundary && v2->dual_intersects_boundary)) {

							double stretch = it->edge_stretch[i]; 

							if (!isolate_edges[gmark].first.is_valid()
								|| (abs(stretch) > abs(isolate_edges[gmark].second)))
							{
								isolate_edges[gmark] = std::make_pair(triEdge(it, i), stretch);
							} 
					}
				}
			}
		}
		
		//insert and remove
		begin_insert();
		for (int i=0; i<isolate_edges.size(); ++i)
		{
			EdgeValue& ev = isolate_edges[i];
			triEdge& e = ev.first;
			Delaunay::Vertex_handle v1 = e.start(); 
			Delaunay::Vertex_handle v2 = e.end();
			vec2 p1 = to_geex(v1->point());
			vec2 p2 = to_geex(v2->point());

			if (ev.second > 0)
			{
				insert((p1+p2)/2);
			} 
			else
			{
				insert((p1+p2)/2);
				remove(p1);
				remove(p2);
			}
		}
		end_insert(false);
	}



	//edit selected

	void Delaunay::split_selected_edge() {
		std::set<Delaunay::Vertex_handle> edgeStart;
		std::vector<vec2> to_insert;
		FOR_EACH_FACE(Delaunay, this, it) {
			if (is_infinite(it)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = (i+1) % 3;
				int j2 = (i+2) % 3;
				Delaunay::Vertex_handle v1 = it->vertex(j1);
				Delaunay::Vertex_handle v2 = it->vertex(j2);
				//if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue;  //exclude boundary vertex

				if (it->selected_edge[i])  //selected edge condition
				{
					if (edgeStart.find(v2) == edgeStart.end())
					{
						edgeStart.insert(v1);
						double xnew = 0.5 * (v1->point().x() + v2->point().x());
						double ynew = 0.5 * (v1->point().y() + v2->point().y());
						to_insert.push_back(vec2(xnew, ynew));
					}
				}
			}
		}
		//insert
		begin_insert() ;
		for (int i=0; i<to_insert.size(); ++i)
		{
			insert(to_insert[i]);
		}
		end_insert(false) ;
		std::cout << to_insert.size() << " points inserted." << std::endl;
	}

	void Delaunay::collapse_selected_edge() {
		std::set<Delaunay::Vertex_handle> to_remove;
		std::vector<vec2> to_insert;
		FOR_EACH_FACE(Delaunay, this, it) {
			if (is_infinite(it)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = (i+1) % 3;
				int j2 = (i+2) % 3;
				Delaunay::Vertex_handle v1 = it->vertex(j1);
				Delaunay::Vertex_handle v2 = it->vertex(j2);
				//if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue;  //exclude boundary vertex

				if (it->selected_edge[i])  //selected edge condition
				{
					if (to_remove.find(v1) == to_remove.end() || to_remove.find(v2) == to_remove.end())
					{
						to_remove.insert(v1);
						to_remove.insert(v2);
						double xnew = 0.5 * (v1->point().x() + v2->point().x());
						double ynew = 0.5 * (v1->point().y() + v2->point().y());
						to_insert.push_back(vec2(xnew, ynew));
					}
				}
			}
		}
		//insert and remove
		begin_insert() ;
		for (int i=0; i<to_insert.size(); ++i)
		{
			insert(to_insert[i]);
		}
		for (auto it=to_remove.begin(); it != to_remove.end(); ++it)
		{
			remove(vec2((*it)->point().x(), (*it)->point().y()));
		}
		end_insert(false) ;
		std::cout << to_remove.size() - to_insert.size() << " points removed." << std::endl;
	}

	void Delaunay::flip_selected_edge() {
		unsigned int count = 0;
		std::vector<triEdge> flipped_edge;

		FOR_EACH_FACE(Delaunay, this, it) {
			for (int i=0; i<3; ++i)
			{
				if (it->selected_edge[i])
				{
					it->selected_edge[i] = false;
					it->neighbor(i)->selected_edge[mirror_index(it, i)] = false;
					flipped_edge.push_back(triEdge(it, (i+1)%3));
					baseclass::flip(it, i);
					count ++;
				}
			}
		}
		//recover flipped edges
		for (int i=0; i<flipped_edge.size(); ++i)
		{
			Delaunay::Face_handle f = flipped_edge[i].face();
			int e = flipped_edge[i].edge_index();
			f->selected_edge[e] = true;
			f->neighbor(e)->selected_edge[mirror_index(f, e)] = true;
		}
		//
		std::cout << count << " edge flipped." << std::endl;
	}

	void Delaunay::update_area(bool info) {
		total_area_ = 0;
		FOR_EACH_VERTEX(Delaunay, this, v) {
			total_area_ += v->area;
		}
		average_area_ = total_area_ / nb_vertices();
		if (info)
		{
			std::cout << "total area = " << total_area_ << std::endl;
			std::cout << "num vertices = " << nb_vertices()  << std::endl;
			std::cout << "average area = " << average_area_ << std::endl;
		}
	}

	//select functions
	void Delaunay::clear_select() {
		FOR_EACH_FACE(Delaunay, this, it) {
			it->selected = false;
			if (this->is_infinite(it)) continue;

			for (int i=0; i<3; ++i)
			{
				it->vertex(i)->selected = false;
				it->selected_edge[i] = false;
			}
		}
	}

	void Delaunay::update_select(const vec2& p) {
		//gx_assert(opened_) ;
		if(Numeric::is_nan(p.x) || Numeric::is_nan(p.y)) {
			std::cerr << "Nan !" << std::endl ;
			return;
		}

		Delaunay::Face_handle f = this->locate(to_cgal(p));
		if (this->is_infinite(f)) return;

		//Face
		if (select_mode_ == FACE) {
			f->selected = !f->selected;
			return;
		}

		if (select_mode_ == EDGE) {
			double min_d = 1e30 ;
			int e = -1;
			for(unsigned int i=0; i<3; i++) {
				Delaunay::Vertex_handle v1 = f->vertex((i+1) % 3);
				double x1 = v1->point().x();
				double y1 = v1->point().y();
				Delaunay::Vertex_handle v2 = f->vertex((i+2) % 3);
				double x2 = v2->point().x();
				double y2 = v2->point().y();
				double len12 = (vec2(x2-x1, y2-y1)).length();
				double cur_d = abs(x2*p.y - y2*p.x - x2*y1 + x1*y2 - x1*p.y + y1*p.x) / len12;
				if (cur_d < min_d)
				{
					min_d = cur_d;
					e = i;
				}
			}

			f->selected_edge[e] = !f->selected_edge[e];
			f->neighbor(e)->selected_edge[this->mirror_index(f, e)] = f->selected_edge[e];
			return;
		}


		if(select_mode_ == SEPARATRIX_SEGMENT) { //user should ensure the initial selected edge is on separatrix and in the right direction
			//select edge e
			double min_d = 1e30 ;
			int e = -1;
			for(unsigned int i=0; i<3; i++) {
				Delaunay::Vertex_handle v1 = f->vertex((i+1) % 3);
				double x1 = v1->point().x();
				double y1 = v1->point().y();
				Delaunay::Vertex_handle v2 = f->vertex((i+2) % 3);
				double x2 = v2->point().x();
				double y2 = v2->point().y();
				double len12 = (vec2(x2-x1, y2-y1)).length();
				double cur_d = abs(x2*p.y - y2*p.x - x2*y1 + x1*y2 - x1*p.y + y1*p.x) / len12;
				if (cur_d < min_d)
				{
					min_d = cur_d;
					e = i;
				}
			}
			//select separatrix segment
			Delaunay::Vertex_handle v = f->vertex((e+1)%3);
			if (v->separatrix_in_degree != 0) //v is on separatrix
			{
				if (f->selected_edge[e])  //unselect separatrix
				{
					Delaunay::Vertex_handle v1 = v;
					Delaunay::Face_handle   f1 = f;
					Delaunay::Vertex_handle v2 = this->incident_vertices(v1, f1);
					//forward
					do 
					{
						for (int i=0; i<3; ++i)
						{
							if (v1 == f1->vertex((i+1)%3))
							{
								f1->selected_edge[i] = false;
								f1->neighbor(i)->selected_edge[this->mirror_index(f1, i)] = false;
							}
						}
						Delaunay::Face_circulator fc = this->incident_faces(v2, f1);
						fc--; fc--;
						f1 = fc;
						v1 = v2;
						v2 = this->incident_vertices(v1, f1);
					} while (v1->singularity == 0 
						&& !v1->dual_intersects_boundary
						&& v1->separatrix_in_degree == 1);
					//backward
					v1 = v;
					f1 = f;
					v2 = this->incident_vertices(v1, f1);
					while (v1->singularity == 0 
						&& !v1->dual_intersects_boundary
						&& v1->separatrix_in_degree == 1)
					{
						Delaunay::Face_circulator fc = this->incident_faces(v1, f1);
						Delaunay::Vertex_circulator vc = this->incident_vertices(v1, f1);
						fc++; fc++;
						vc++; vc++; vc++;
						f1 = fc;
						v2 = v1;
						v1 = vc;
						for (int i=0; i<3; ++i)
						{
							if (v1 == f1->vertex((i+1)%3))
							{
								f1->selected_edge[i] = false;
								f1->neighbor(i)->selected_edge[this->mirror_index(f1, i)] = false;
							}
						}
					}
				}
				else { //select separatrix
					Delaunay::Vertex_handle v1 = v;
					Delaunay::Face_handle   f1 = f;
					Delaunay::Vertex_handle v2 = this->incident_vertices(v1, f1);
					//forward
					do 
					{
						for (int i=0; i<3; ++i)
						{
							if (v1 == f1->vertex((i+1)%3))
							{
								f1->selected_edge[i] = true;
								f1->neighbor(i)->selected_edge[this->mirror_index(f1, i)] = true;
							}
						}
						Delaunay::Face_circulator fc = this->incident_faces(v2, f1);
						fc--; fc--;
						f1 = fc;
						v1 = v2;
						v2 = this->incident_vertices(v1, f1);
					} while (v1->singularity == 0 
						&& !v1->dual_intersects_boundary
						&& v1->separatrix_in_degree == 1);
					//backward
					v1 = v;
					f1 = f;
					v2 = this->incident_vertices(v1, f1);
					while (v1->singularity == 0 
						&& !v1->dual_intersects_boundary
						&& v1->separatrix_in_degree == 1)
					{
						Delaunay::Face_circulator fc = this->incident_faces(v1, f1);
						Delaunay::Vertex_circulator vc = this->incident_vertices(v1, f1);
						fc++; fc++;
						vc++; vc++; vc++;
						f1 = fc;
						v2 = v1;
						v1 = vc;
						for (int i=0; i<3; ++i)
						{
							if (v1 == f1->vertex((i+1)%3))
							{
								f1->selected_edge[i] = true;
								f1->neighbor(i)->selected_edge[this->mirror_index(f1, i)] = true;
							}
						}
					}
				}
			}
		}

		if (select_mode_ == VERTEX)
		{
			double min_d = 1e30 ;
			Delaunay::Vertex_handle v = 0 ;
			for(unsigned int i=0; i<3; i++) {
				Delaunay::Vertex_handle v1 = f->vertex(i);
				double cur_d = (to_geex(v1->point()) - p).length2() ;
				if(cur_d < min_d) {
					min_d = cur_d ;
					v = v1 ;
				}

			}
			v->selected = !v->selected;
			return;
		}

	}


	///dxy: about boundary clip ... rarely used functions
	Polygon_2 Delaunay::to_cgal_polygon(const Polygon2& png)
	{
		std::vector<Point> pVec(png.size());
		for (int i=0; i<png.size(); ++i)
		{
			const PolygonEdge& e = png[i];
			const PolygonVertex& v = e.vertex[0];
			pVec[i] = to_cgal(v);
		}
		Polygon_2 to(pVec.begin(), pVec.end());
		return to;
	}

	double Delaunay::boundary_clipped_length(const CGAL::Object& obj)
	{
		double length = 0;
		//input is ray
		if (const Ray_2 *ray = CGAL::object_cast<Ray_2>(&obj))
		{
			vec2 p1 = to_geex(ray->source());
			if (in_boundary(p1))
			{
				length = std::numeric_limits<double>::max();
				bool has_intersect = false;  ///test
				for (int i=0; i<cgal_boundary_.size(); ++i)
				{
					Segment_2 edge = cgal_boundary_.edge(i);
					CGAL::Object intersect = CGAL::intersection(*ray, edge);
					if (const Point* point = CGAL::object_cast<Point>(&intersect))
					{
						has_intersect = true; ///test
						vec2 p2 = to_geex(*point); 
						vec2 p1p2 = p2 - p1;
						length = min(length, p1p2.length());
					}
					else if (const Segment_2* segment = CGAL::object_cast<Segment_2>(&intersect))
					{
						has_intersect = true; ///test
						vec2 p2 = to_geex(segment->source());
						vec2 p1p2 = p2 - p1;
						double len12 = p1p2.length();
						vec2 p3 = to_geex(segment->target());
						vec2 p1p3 = p3 - p1;
						double len13 = p1p3.length();
						double len = min(len12, len13);
						length = min(len, length);
					}
				}
				///test
				if (!has_intersect && length > 10)
				{
					return 12.1;
				}
				///
			}
		}
		//input is segment
		else if (const Segment_2 *seg = CGAL::object_cast<Segment_2>(&obj))
		{
			vec2 p1 = to_geex(seg->source());
			vec2 p2 = to_geex(seg->target());
			bool in1 = in_boundary(p1);
			bool in2 = in_boundary(p2);
			if (in1 && in2)  //assume: convex boundary
			{
				vec2 p1p2 = p2 - p1;
				return p1p2.length();

			}
			if (!in1 && !in2)  
			{
				return 0;
			}
			if (in2)
			{
				p1 = p2;
			}
			//
			length = std::numeric_limits<double>::max();
			for (int i=0; i<cgal_boundary_.size(); ++i)
			{
				Segment_2 edge = cgal_boundary_.edge(i);
				CGAL::Object intersect = CGAL::intersection(*seg, edge);
				if (const Point* point = CGAL::object_cast<Point>(&intersect))
				{
					p2 = to_geex(*point); 
					vec2 p1p2 = p2 - p1;
					length = min(length, p1p2.length());
				}
				else if (const Segment_2* segment = CGAL::object_cast<Segment_2>(&intersect))
				{
					p2 = to_geex(segment->source());
					vec2 p1p2 = p2 - p1;
					double len12 = p1p2.length();
					vec2 p3 = to_geex(segment->target());
					vec2 p1p3 = p3 - p1;
					double len13 = p1p3.length();
					double len = min(len12, len13);
					length = min(len, length);
				}
			}
		}
		return length;
	}

	///class triEdge
	
	///dxy add end


}
