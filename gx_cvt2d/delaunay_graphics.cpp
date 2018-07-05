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

#include "delaunay_graphics.h"
#include <glut_viewer/glut_viewer.h>
#include <delaunay_cvt.h>

namespace Geex {


	static const double c1 = 0.35 ;
	static const double c2 = 0.5 ;
	static const double c3 = 1.0 ;

	static double color_table[12][3] = 
	{
		{c3, c2, c2},
		{c2, c3, c2},
		{c2, c2, c3},
		{c2, c3, c3},
		{c3, c2, c3},
		{c3, c3, c2},

		{c1, c2, c2},
		{c2, c1, c2},
		{c2, c2, c1},
		{c2, c1, c1},
		{c1, c2, c1},
		{c1, c1, c2}

	} ;

	static int random_color_index_ = 0 ;


	static void gl_random_color() {
		glColor3f(
			color_table[random_color_index_][0], 
			color_table[random_color_index_][1], 
			color_table[random_color_index_][2]
		) ;
		random_color_index_ = (random_color_index_ + 1) % 12 ;
	}

	static void gl_random_color(int index) {
		random_color_index_ = index % 12 ;
		gl_random_color() ;
	}

	static void gl_randomize_colors() {
		random_color_index_ = 0 ;
	}

	///dxy copy: from 3D CVT
	static void gl_vertex_color(int degree) {
		if(degree==6) {
			glColor3f(213/255.0, 222/255.0, 157/255.0) ;
		}
		else if(degree==7) {
			glColor3f(125/255.0, 176/255.0, 223/255.0) ;
		}
		else if(degree==5) {
			glColor3f(192/255.0, 121/255.0, 165/255.0) ;
		}
		else if(degree>7) {
			glColor3f(98/255.0, 57/255.0, 115/255.0) ;
		}
		else {
			glColor3f(0, 0.25, 0.5) ;
			glColor3f(41/255.0, 118/255.0, 102/255.0) ;
		}
	}
	///dxy add
	//hsv ([0,360], [0,1], [0,1]),  rgb ([0,1]...)
	static void HSVtoRGB(real& fR, real& fG, real& fB, real& fH, real& fS, real& fV) {
		real fC = fV * fS; // Chroma
		real fHPrime = fmod(fH / 60.0, 6);
		real fX = fC * (1 - fabs(fmod(fHPrime, 2) - 1));
		real fM = fV - fC;

		if(0 <= fHPrime && fHPrime < 1) {
			fR = fC;
			fG = fX;
			fB = 0;
		} else if(1 <= fHPrime && fHPrime < 2) {
			fR = fX;
			fG = fC;
			fB = 0;
		} else if(2 <= fHPrime && fHPrime < 3) {
			fR = 0;
			fG = fC;
			fB = fX;
		} else if(3 <= fHPrime && fHPrime < 4) {
			fR = 0;
			fG = fX;
			fB = fC;
		} else if(4 <= fHPrime && fHPrime < 5) {
			fR = fX;
			fG = 0;
			fB = fC;
		} else if(5 <= fHPrime && fHPrime < 6) {
			fR = fC;
			fG = 0;
			fB = fX;
		} else {
			fR = 0;
			fG = 0;
			fB = 0;
		}

		fR += fM;
		fG += fM;
		fB += fM;
	}
	static void hsv2rgb(vec3& hsv, vec3& rgb) 
	{
		real& fR = rgb[0];
		real& fG = rgb[1];
		real& fB = rgb[2];
		real& fH = hsv[0];
		real& fS = hsv[1];
		real& fV = hsv[2];
		HSVtoRGB(fR, fG, fB, fH, fS, fV);
	}
	///dxy add end

	DelaunayGraphics::DelaunayGraphics(DelaunayCVT* CVT) : delaunay_(CVT->delaunay()), CVT_(CVT) {
		show_domain_ = GL_TRUE ; 
		vertices_size_ = 0.0 ;
		centers_size_ = 0.0 ;
		show_primal_mesh_ = GL_FALSE ;
		show_dual_mesh_ = GL_TRUE ;
		colorize_ = GL_TRUE ;
		show_primal_ = GL_FALSE ;
		show_cells_ = GL_FALSE ;
		show_energy_ = GL_FALSE ;
		show_field_ = GL_FALSE ;
		quad_ratio_ = 0.0 ;
		show_affinity_ = GL_FALSE ;

		///dxy add
		show_lloyd_grad_ = GL_FALSE;
		show_direction_grad_ = GL_FALSE;
		show_constrained_edge_ = GL_FALSE;
		lloyd_grad_magnification_ = 1;
		direct_grad_magnification_ = 1;
		total_grad_magnification_ = 1;
		stretch_long_rank_  = 0;
		stretch_short_rank_ = 0;
		stretch_all_rank_ = -1;
		///
	}

	void DelaunayGraphics::draw() {
		non_convex_ = delaunay_->non_convex_mode_ ;
		glDisable(GL_LIGHTING) ;
		glUseProgramObjectARB(0) ;

		gl_randomize_colors() ;

		if(show_domain_) {
			draw_domain() ;
		}

		if(show_affinity_) {
			draw_affinity() ;
		}


		///dxy add
		if(show_cells_) {
			draw_non_hex_cells() ;
		}

		if(show_energy_) {
			draw_energy() ;
		}

		if (show_relative_area_)
		{
			draw_relative_area();
		}
		if (show_triangle_area_)
		{
			draw_primal_area();
		}
		///


		if(show_field_) {
			draw_field() ;
		}



		glDisable(GL_TEXTURE_2D) ;
		glDisable(GL_TEXTURE_1D) ;
		glUseProgramObjectARB(0) ;

		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST) ;
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
		glEnable(GL_LINE_SMOOTH) ;
		glEnable(GL_POINT_SMOOTH) ;
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
		glEnable(GL_BLEND) ;
		glDepthMask(GL_FALSE) ;


		if(show_primal_mesh_) {
			draw_primal_mesh() ;
		}

		if(show_dual_mesh_) {
			draw_dual_mesh() ;
		}

		if(centers_size_ != 0.0) {
			draw_centers() ;
		}

		if(vertices_size_ != 0.0) {
			draw_vertices() ;
		}

		///dxy add
		if (show_constrained_edge_)
		{
			draw_constrained_edge();
		}

		if (show_lloyd_grad_)
		{
			draw_lloyd_grad();
		}

		if (show_direction_grad_)
		{
			draw_direction_grad();
		}

		if (show_total_grad_)
		{
			draw_total_grad();
		}

		if (show_direction_field_)
		{
			draw_direction_field();
		}

		if (show_long_edge_)
		{
			draw_isolate_long_edge();
		}
		if (show_short_edge_)
		{
			draw_isolate_short_edge();
		}

		if (show_edge_stretch_)
		{
			draw_edge_stretch();
		}

		if (show_separatrix_)
		{
			draw_separatrix();
		}

		if (show_selected_)
		{
			draw_selected();
		}
		///



		glDisable(GL_LINE_SMOOTH) ;
		glDisable(GL_BLEND) ;
		glDepthMask(GL_TRUE) ;
	}

	void DelaunayGraphics::draw_domain() {
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glBegin(GL_LINES) ;
		for(unsigned int i=0; i<delaunay_->boundary_.size(); i++) {
			glVertex(delaunay_->boundary_[i].vertex[0]) ;
			glVertex(delaunay_->boundary_[i].vertex[1]) ;
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_affinity() {
		if(CVT_->has_affinity()) {
			glLineWidth(1) ;
			glColor3f(0.0, 0.0, 0.5) ;
			glBegin(GL_LINES) ;
			for(unsigned int i=0; i<256; i++) {
				for(unsigned int j=i+1; j<256; j++) {
					if(CVT_->affinity(i,j) < 0.0) {
						glVertex(delaunay_->vertex(i)) ;
						glVertex(delaunay_->vertex(j)) ;
					}
				}
			}
			glEnd() ;
		}
	}

	void DelaunayGraphics::draw_vertices() {
		glDisable(GL_LIGHTING) ;
		glPointSize(int(vertices_size_ * 20)) ;
		int w,h ;
		glut_viewer_get_screen_size(&w, &h) ;
		notify_screen_resize(w,h) ;
		glEnable(GL_POINT_SMOOTH) ;
		glBegin(GL_POINTS) ;
		//        begin_spheres() ;
		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			if(v->locked) {
				glColor3f(1.0,0.0,0.0) ;
			} else {
				glColor3f(0.0,0.0,0.0) ;
			}
			glVertex2f(v->point().x(), v->point().y()) ;
		}
		//        end_spheres() ;
		glEnd() ;
		if(vertices_size_ > 0.6) {
			glBegin(GL_LINES) ;
			FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
				double R = 0.5 * radius(v) ;
				vec2 p = to_geex(v->point()) ;

				vec2 U,V ;
				CVT_->query_anisotropy(p,U,V) ;
				/*
				vec2 U(cos(v->theta), sin(v->theta)) ;
				U *= R ;
				vec2 V(-U.y, U.x) ; 
				*/              

				U *= 0.3 ;
				V *= 0.3 ;
				glVertex(p-U) ; glVertex(p+U) ;
				glVertex(p-V) ; glVertex(p+V) ;
			}
			glEnd() ;
		}
	}

	void DelaunayGraphics::draw_centers() {
		glDisable(GL_LIGHTING) ;
		glPointSize(int(centers_size_ * 20)) ;
		glColor3f(0.0,0.5,0.0) ;
		int w,h ;
		glut_viewer_get_screen_size(&w, &h) ;
		notify_screen_resize(w,h) ;
		glEnable(GL_POINT_SMOOTH) ;
		glBegin(GL_POINTS) ;
		//        begin_spheres() ;
		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			double V ;
			vec2 g ;
			CVT_->get_cell_centroid(v, g, V) ;
			vec2 p = g ;
			glVertex2f(p.x, p.y) ;
		}
		glEnd() ;
		//        end_spheres() ;
	}

	void DelaunayGraphics::draw_cells(bool colorize, bool mesh) {
		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			if(colorize) {gl_random_color() ; }
			if(delaunay_->dual_cell_intersects_boundary(v)) { 
				Polygon2* P = delaunay_->dual_convex_clip(v) ;
				if(mesh) {
					glBegin(GL_LINES) ;
					for(unsigned int i=0; i<P->size(); i++) {
						glVertex((*P)[i].vertex[0]) ;
						glVertex((*P)[i].vertex[1]) ;
					}
					glEnd() ;
				} else {
					glBegin(GL_TRIANGLES) ;
					for(unsigned int i=0; i<P->size(); i++) {
						glVertex(to_geex(v->point())) ;
						glVertex((*P)[i].vertex[0]) ;
						glVertex((*P)[i].vertex[1]) ;
					}
					glEnd() ;
				}
			} else {
				glBegin(GL_POLYGON) ;
				Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
				do {
					glVertex(it->dual) ;
					it++ ;
				} while(it != delaunay_->incident_faces(v)) ;
				glEnd() ;
			} 
		}
	}

	//////////////////////////////////////////////////////////////////////////
	///dxy change: add info
	//void DelaunayGraphics::draw_non_hex_cells() {
	void DelaunayGraphics::draw_non_hex_cells(bool info) {
		unsigned nlt4 = 0;
		unsigned n4 = 0;
		unsigned n5 = 0;
		unsigned n7 = 0;
		unsigned n8 = 0;
		unsigned ngt8 = 0;
		unsigned nsing = 0; //num singularity
		unsigned nwsing = 0;  //num weighted singularity
		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			if(delaunay_->dual_cell_intersects_boundary(v)) { 
				//Polygon2* P = delaunay_->dual_convex_clip(v) ;
				/////
				//glColor4f(0.0, 0.5, 0.0, alpha);
				////glColor3f(0.0, 0.5, 0.0) ;
				//glBegin(GL_TRIANGLES) ;
				//for(unsigned int i=0; i<P->size(); i++) {
				//	glVertex(to_geex(v->point())) ;
				//	glVertex((*P)[i].vertex[0]) ;
				//	glVertex((*P)[i].vertex[1]) ;
				//}
				//glEnd() ;
			} else  {
				int valence = delaunay_->dual_facet_degree(v);
				if(valence != 6) {
					nsing++;
					nwsing += abs(valence - 6);
					if (valence < 4)
					{
						nlt4++;
						break;
					}
					if (valence > 8)
					{
						ngt8++;
						break;
					}
					switch (valence)
					{
					case 4:
						n4++;
						break;
					case 5:
						n5++;
						break;
					case 7:
						n7++;
						break;
					case 8:
						n8++;
						break;
					default:
						break;
					}
					///dxy change
					//glColor3f(0.0, 0.0, 0.5) ;
					gl_vertex_color(valence);
					///
					glBegin(GL_POLYGON) ;
					Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
					do {
						glVertex(it->dual) ;
						it++ ;
					} while(it != delaunay_->incident_faces(v)) ;
					glEnd() ;
				}
			} 
		}

		if (info)
		{
			std::cout << "nlt4 = " << nlt4 << ", n4 = " << n4 << std::endl;
			std::cout << "n5 = " << n5 << ", n7 = " << n7 << std::endl;
			std::cout << "n8 = " << n8 << ", ngt8 = " << ngt8 << std::endl;
			std::cout << "ns = " << nsing << std::endl;
			std::cout << "nws = " << nwsing << std::endl; 
		}
	}
	////////////////////////////////////////////////////////////////////////// dxy change end

	inline bool is_quad(Delaunay* del, Delaunay::Face_handle f, int i, double ratio) {
		Delaunay::Face_handle g = f->neighbor(i) ;
		if(del->is_infinite(g)) { return false ; }

		vec2 p1 = to_geex(f->vertex(0)->point()) ;
		vec2 q1 = to_geex(g->vertex(0)->point()) ;

		vec2 c1 = f->dual ;
		vec2 c2 = g->dual ;

		double R1 = (p1 - c1).length() ;
		double R2 = (q1 - c2).length() ;
		double R = 0.5 * (R1 + R2) ;

		double dC = (c2-c1).length() ;
		return (dC < ratio * R) ;
	}

	void DelaunayGraphics::draw_primal_mesh() {
		glLineWidth(1) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_LINES) ;
		FOR_EACH_FACE(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) { continue ; }

			for(unsigned int i=0; i<3; i++) {
				unsigned int j1 = i + 1 ; 
				if(j1 == 3) { j1 = 0 ; }
				unsigned int j2 = (j1 + 1) ;
				if(j2 == 3) { j2 = 0 ; }

				if(!is_quad(delaunay_, it, i, quad_ratio_)) {
					Delaunay::Vertex_handle v1 = it->vertex(j1) ;
					Delaunay::Vertex_handle v2 = it->vertex(j2) ;
					glVertex2f(v1->point().x(), v1->point().y()) ;
					glVertex2f(v2->point().x(), v2->point().y()) ;
				}
			}
		}
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
	}

	void DelaunayGraphics::draw_dual_mesh() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glLineWidth(9) ;
		glColor3f(1.0, 1.0, 1.0) ;
		draw_cells(false, true) ;
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 0.0) ;
		draw_cells(false, true) ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
	}

	void DelaunayGraphics::is_min_max(Delaunay::Vertex_handle v, bool& is_min, bool& is_max) {
		is_min = true ;
		is_max = true ;
		Delaunay::Vertex_circulator it = delaunay_->incident_vertices(v) ;
		do {
			if(!delaunay_->is_infinite(it)) {
				is_min = is_min && v->energy < it->energy ;
				is_max = is_max && v->energy > it->energy ;
			}
			it++ ;
		} while(it != delaunay_->incident_vertices(v)) ;
	}

	//////////////////////////////////////////////////////////////////////////
	///dxy change
	void DelaunayGraphics::draw_energy() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		double emin = 1e30 ;
		double emax = -1e30 ;

		Delaunay::Vertex_handle vmin = 0 ;
		Delaunay::Vertex_handle vmax = 0 ;

		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			if(v->energy < emin) { emin = v->energy; vmin = v ; }
			if(v->energy > emax) { emax = v->energy; vmax = v ; }
		}
		//        std::cerr << "emin = " << emin << " emax =" << emax << std::endl ;

		double emean = 0.5 * (emin + emax) ;

		/* if(emax - emean > emean - emin) */ {
			//glColor3f(1.0, 0.5, 0.5) ;
			glColor3f(1.0, 0.0, 0.0);
			if(vmax != 0) { draw_dual_facet(vmax) ; }
		} 
		/* else */ {

			//glColor3f(0.5, 0.5, 1.0) ;
			glColor3f(0.0, 1.0, 0.0);
			if(vmin != 0) { draw_dual_facet(vmin) ; }
		}

		double scale = (emax - emin) ;
		if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
		scale = 1.0 / scale ;

		///dxy change
		const vec3 c1(120, 1.0, 1.0);
		const vec3 c2(0.0, 1.0, 1.0);
		///

		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			double e = scale * (v->energy - emin) ;
			///dxy add
			vec3 mix_hsv = mix(c1, c2, e);
			vec3 mix_rgb;
			hsv2rgb(mix_hsv, mix_rgb);
			glColor3f(mix_rgb.x, mix_rgb.y, mix_rgb.z);
			///
			draw_dual_facet(v) ;
		}
	}
	////////////////////////////////////////////////////////////////////////// dxy change end

	void DelaunayGraphics::draw_dual_facet(Delaunay::Vertex_handle v) {
		if(delaunay_->dual_cell_intersects_boundary(v)) { 
			Polygon2* P = delaunay_->dual_convex_clip(v) ;
			glBegin(GL_TRIANGLES) ;
			for(unsigned int i=0; i<P->size(); i++) {
				glVertex(to_geex(v->point())) ;
				glVertex((*P)[i].vertex[0]) ;
				glVertex((*P)[i].vertex[1]) ;
			}
			glEnd() ;
		} else  {
			glBegin(GL_POLYGON) ;
			Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			do {
				glVertex(it->dual) ;
				it++ ;
			} while(it != delaunay_->incident_faces(v)) ;
			glEnd() ;
		} 
	}


	void DelaunayGraphics::draw_field() {
		double x_min, y_min, z_min ;
		double x_max, y_max, z_max ;
		delaunay_->get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
		double dx = x_max - x_min ;
		double dy = y_max - y_min ;
		glPointSize(4) ;
		glBegin(GL_POINTS) ;



		for(int i=0; i<100; i++) {
			for(int j=0; j<100; j++) {
				vec2 p(x_min + dx * double(i)/100.0, y_min + dy * double(j)/100.0) ;
				gl_random_color(delaunay_->segments().locate(p)) ;
				glVertex(p) ;
			}
		}

		glPointSize(8) ;
		for(SegmentDelaunay::Vertex_iterator it = delaunay_->segments().finite_vertices_begin() ;
			it != delaunay_->segments().finite_vertices_end(); it++
			) {
				gl_random_color(it->index) ;
				glVertex(to_geex(it->point())) ;
		}



		glEnd() ;



	}

	double DelaunayGraphics::radius(Delaunay::Vertex_handle v) {
		double result = 1e30 ;
		vec2 p1 = to_geex(v->point()) ;
		Delaunay::Vertex_circulator it = delaunay_->incident_vertices(v) ;
		do {
			if(!delaunay_->is_infinite(it)) {
				vec2 p2 = to_geex(it->point()) ;
				result = gx_min(result, (p2 - p1).length2()) ;
			}
			it++ ;
		} while(it != delaunay_->incident_vertices(v)) ;
		return sqrt(result) ;
	}



	//////////////////////////////////////////////////////////////////////////
	///dxy add
	void DelaunayGraphics::draw_constrained_edge() {
		glLineWidth(1) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_LINES) ;
		FOR_EACH_FACE(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) { continue ; }

			for(unsigned int i=0; i<3; i++) {
				unsigned int j1 = i + 1 ; 
				if(j1 == 3) { j1 = 0 ; }
				unsigned int j2 = (j1 + 1) ;
				if(j2 == 3) { j2 = 0 ; }

				if(!is_quad(delaunay_, it, i, quad_ratio_)) {
					Delaunay::Vertex_handle v1 = it->vertex(j1) ;
					Delaunay::Vertex_handle v2 = it->vertex(j2) ;
					int cw_j1 = (j1 + 2) % 3;
					Delaunay::Face_handle f = it->neighbor(cw_j1);
					if (delaunay_->is_infinite(f))
					{
						glColor3f(0.0, 0.0, 0.0);   //boundary edge: black
						glVertex2f(v1->point().x(), v1->point().y()) ;
						glVertex2f(v2->point().x(), v2->point().y()) ;
					}
					else {
						vec2 p1 = it->dual;
						vec2 p2 = f->dual;
						if (!delaunay_->in_boundary(p1) || !delaunay_->in_boundary(p2))
						{
							glColor3f(1.0, 0.0, 0); //constrained edge: red
							glVertex2f(v1->point().x(), v1->point().y()) ;
							glVertex2f(v2->point().x(), v2->point().y()) ;
						}
					}
				}
			}
		}
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
	}

	void DelaunayGraphics::draw_lloyd_grad() {
		glLineWidth(3) ;
		glColor3f(1.0, 0.0, 0.0) ;
		FOR_EACH_VERTEX(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) {continue ; }
			vec2 gstart = to_geex(it->point()) ;
			vec2 grad = -1/*10000*/ * lloyd_grad_magnification_ * it->lloyd_grad;
			draw_vector(grad, gstart);
		}
	}

	void DelaunayGraphics::draw_direction_grad() {
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 1.0) ;
		glBegin(GL_LINES) ;
		FOR_EACH_VERTEX(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) {continue ; }
			vec2 gstart = to_geex(it->point()) ;
			vec2 grad = -1/*10000*/ * direct_grad_magnification_ * it->direction_grad ;
			draw_vector(grad, gstart);
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_total_grad() {
		glLineWidth(3) ;
		glColor3f(0.0, 1.0, 0.0) ;
		glBegin(GL_LINES) ;
		FOR_EACH_VERTEX(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) {continue ; }
			vec2 gstart = to_geex(it->point()) ;
			vec2 grad = -1/*10000*/ * total_grad_magnification_ * (it->direction_grad + it->lloyd_grad) ;
			draw_vector(grad, gstart);
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_direction_field() {
		double average_length = sqrt(delaunay_->average_area_);
		glLineWidth(2.5);
		glColor3f(0.0, 0.0, 0.0);
		FOR_EACH_VERTEX(Delaunay, delaunay_, it) {
			vec2 start = to_geex(it->point());
			double field_angle;
			vec2 field_grad;
			CVT_->get_field(start, field_angle, field_grad);
			vec2 field(cos(field_angle), sin(field_angle));
			field *= 0.7 * average_length;
			draw_vector(field, start);
		}
	}

	//
	void DelaunayGraphics::draw_relative_area() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		double emin = 1e30 ;
		double emax = -1e30 ;

		Delaunay::Vertex_handle vmin = 0 ;
		Delaunay::Vertex_handle vmax = 0 ;

		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			if(v->area < emin) { emin = v->area; vmin = v ; }
			if(v->area > emax) { emax = v->area; vmax = v ; }
		}

		double emean = 0.5 * (emin + emax) ;

		/* if(emax - emean > emean - emin) */ {
			glColor3f(1.0, 0.0, 1.0);
			if(vmax != 0) { draw_dual_facet(vmax) ; }
		} 
		/* else */ {
			glColor3f(0.0, 1.0, 1.0);
			if(vmin != 0) { draw_dual_facet(vmin) ; }
		}

		double scale = (emax - emin) ;
		if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
		scale = 1.0 / scale ;

		const vec3 c1(120, 1.0, 1.0);
		const vec3 c2(0.0, 1.0, 1.0);

		FOR_EACH_VERTEX(Delaunay, delaunay_, v) {
			double e = scale * (v->area - emin) ;
			vec3 mix_hsv = mix(c1, c2, e);
			vec3 mix_rgb;
			hsv2rgb(mix_hsv, mix_rgb);
			glColor3f(mix_rgb.x, mix_rgb.y, mix_rgb.z);
			draw_dual_facet(v) ;
		}
	}

	void DelaunayGraphics::draw_primal_area() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		double emin = 1e30 ;
		double emax = -1e30 ;

		Delaunay::Face_handle fmin = 0 ;
		Delaunay::Face_handle fmax = 0 ;

		FOR_EACH_FACE(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;

			double farea = area(f->vertex(0)->point(),
				f->vertex(1)->point(),
				f->vertex(2)->point()
				);
			if(farea < emin) { emin = farea; fmin = f ; }
			if(farea > emax) { emax = farea; fmax = f ; }
		}
		double emean = 0.5 * (emin + emax) ;

		/* if(emax - emean > emean - emin) */ {
			glColor3f(1.0, 0.0, 1.0);
			if(fmax != 0) { draw_primal_triangle(fmax) ; }
		} 
		/* else */ {
			glColor3f(0.0, 1.0, 1.0);
			if(fmin != 0) { draw_primal_triangle(fmin) ; }
		}

		double scale = (emax - emin) ;
		if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
		scale = 1.0 / scale ;

		const vec3 c1(120, 1.0, 1.0);
		const vec3 c2(0.0, 1.0, 1.0);

		FOR_EACH_FACE(Delaunay, delaunay_, f) {
			double farea = area(f->vertex(0)->point(),
				f->vertex(1)->point(),
				f->vertex(2)->point()
				);
			double e = scale * (farea - emin) ;
			vec3 mix_hsv = mix(c1, c2, e);
			vec3 mix_rgb;
			hsv2rgb(mix_hsv, mix_rgb);
			glColor3f(mix_rgb.x, mix_rgb.y, mix_rgb.z);
			draw_primal_triangle(f) ;
		}
	}

	//
	void DelaunayGraphics::draw_long_edge()
	{
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = i % 3;
				int j2 = (j1 + 1) % 3;
				Delaunay::Vertex_handle v1 = f->vertex(j1);
				Delaunay::Vertex_handle v2 = f->vertex(j2);
				if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue;  //exclude boundary vertex
				int d1 = delaunay_->dual_facet_degree(v1);
				int d2 = delaunay_->dual_facet_degree(v2);

				if (d1 + d2 >= 14) //long edge condition
				{
					count++;
					glVertex2f(v1->point().x(), v1->point().y());
					glVertex2f(v2->point().x(), v2->point().y());
				}
			}
		}
		glEnd();
		std::cout << count/2 << " long edges." << std::endl;
	}

	void DelaunayGraphics::draw_short_edge()
	{
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.5, 0.0, 0.0) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = i % 3;
				int j2 = (j1 + 1) % 3;
				Delaunay::Vertex_handle v1 = f->vertex(j1);
				Delaunay::Vertex_handle v2 = f->vertex(j2);
				if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue; //exclude boundary vertex
				int d1 = delaunay_->dual_facet_degree(v1);
				int d2 = delaunay_->dual_facet_degree(v2);

				if (d1 + d2 <= 10) //short edge condition
				{
					count ++;
					glVertex2f(v1->point().x(), v1->point().y());
					glVertex2f(v2->point().x(), v2->point().y());
				}
			}
		}
		glEnd();
		std::cout << count/2 << " short edges." << std::endl;
	}

	void DelaunayGraphics::draw_isolate_long_edge()
	{
		std::set<Delaunay::Vertex_handle> incident_vertices; //vertices incident to a long edge
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = i % 3;
				int j2 = (j1 + 1) % 3;
				Delaunay::Vertex_handle v1 = f->vertex(j1);
				Delaunay::Vertex_handle v2 = f->vertex(j2);
				if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue;  //exclude boundary vertex
				int d1 = delaunay_->dual_facet_degree(v1);
				int d2 = delaunay_->dual_facet_degree(v2);

				if (d1 + d2 >= 14) //long edge condition
				{
					if (incident_vertices.find(v1) != incident_vertices.end() ||
						incident_vertices.find(v2) != incident_vertices.end()) continue; //isolate edge condition
					count++;
					glVertex2f(v1->point().x(), v1->point().y());
					glVertex2f(v2->point().x(), v2->point().y());
					incident_vertices.insert(v1);
					incident_vertices.insert(v2);
				}
			}
		}
		glEnd();
		std::cout << count << " long edges." << std::endl;
	}

	void DelaunayGraphics::draw_isolate_short_edge()
	{
		std::set<Delaunay::Vertex_handle> incident_vertices; //vertices incident to a short edge
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.5, 0.0, 0.0) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				int j1 = i % 3;
				int j2 = (j1 + 1) % 3;
				Delaunay::Vertex_handle v1 = f->vertex(j1);
				Delaunay::Vertex_handle v2 = f->vertex(j2);
				if (v1->dual_intersects_boundary || v2->dual_intersects_boundary) continue; //exclude boundary vertex
				int d1 = delaunay_->dual_facet_degree(v1);
				int d2 = delaunay_->dual_facet_degree(v2);

				if (d1 + d2 <= 10) //short edge condition
				{
					if (incident_vertices.find(v1) != incident_vertices.end() ||
						incident_vertices.find(v2) != incident_vertices.end()) continue; //isolate edge condition
					count ++;
					glVertex2f(v1->point().x(), v1->point().y());
					glVertex2f(v2->point().x(), v2->point().y());
					incident_vertices.insert(v1);
					incident_vertices.insert(v2);
				}
			}
		}
		glEnd();
		std::cout << count << " short edges." << std::endl;

	}

	void DelaunayGraphics::draw_separatrix() {
		delaunay_->calc_vertex_singularity();

		//collect irregular vertices && clear separatrix mark
		std::vector<Delaunay::Vertex_handle> irregular_vertices;
		FOR_EACH_VERTEX(Delaunay, delaunay_, it) {
			it->separatrix_in_degree = 0;
			if (it->singularity != 0)
			{
				irregular_vertices.push_back(it);
			}
		}

		//collect separatrix edges
		typedef std::pair<Delaunay::Vertex_handle, Delaunay::Face_handle> dirEdge;  //<v,f> : Edge(v, incident_vertices(v,f))

		std::vector<dirEdge> separatrix0;
		for (int i=0; i<irregular_vertices.size(); ++i)
		{
			Delaunay::Vertex_handle v = irregular_vertices[i];
			int degree = delaunay_->dual_facet_degree(v);
			Delaunay::Face_circulator f = delaunay_->incident_faces(v);
			std::queue<dirEdge> Q;
			do 
			{
				Q.push(std::make_pair(v, f));
				f++;
			} while (f != delaunay_->incident_faces(v));

			while (!Q.empty())
			{
				dirEdge top = Q.front();
				Q.pop();
				separatrix0.push_back(top);

				Delaunay::Vertex_handle v1 = top.first;
				Delaunay::Face_handle   f1 = top.second;
				Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, f1);
				v2->separatrix_in_degree += 1;

				if (!v2->dual_intersects_boundary 
					&& v2->singularity == 0)
					//&& v2->separatrix_in_degree < 2)
				{
					Delaunay::Face_circulator fc = delaunay_->incident_faces(v2, f1);
					fc--; fc--;
					Q.push(std::make_pair(v2, fc));
				}
			}
		}

		//
		//std::vector<dirEdge> separatrix1;
		//for (int i=0; i<irregular_vertices.size(); ++i)
		//{
		//	Delaunay::Vertex_handle v = irregular_vertices[i];
		//	int degree = delaunay_->dual_facet_degree(v);
		//	Delaunay::Face_circulator f = delaunay_->incident_faces(v);
		//	std::queue<dirEdge> Q;
		//	do 
		//	{
		//		Q.push(std::make_pair(v, f));
		//		f++;
		//	} while (f != delaunay_->incident_faces(v));

		//	while (!Q.empty())
		//	{
		//		dirEdge top = Q.front();
		//		Q.pop();
		//		separatrix1.push_back(top);

		//		Delaunay::Vertex_handle v1 = top.first;
		//		Delaunay::Face_handle   f1 = top.second;
		//		Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, f1);
		//		//v2->separatrix_in_degree += 1;

		//		if (!v2->dual_intersects_boundary 
		//			&& v2->singularity == 0
		//			&& v2->separatrix_in_degree < 2)
		//		{
		//			Delaunay::Face_circulator fc = delaunay_->incident_faces(v2, f1);
		//			fc--; fc--;
		//			Q.push(std::make_pair(v2, fc));
		//		}
		//	}
		//}

		//draw separatrix
		glLineWidth(3);
		glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		for(int i=0; i<separatrix0.size(); ++i) {
			dirEdge& edge = separatrix0[i];
			Delaunay::Vertex_handle v1 = edge.first;
			Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, edge.second);
			glVertex2f(v1->point().x(), v1->point().y());
			glVertex2f(v2->point().x(), v2->point().y());
		}
		//
		/*glColor3f(0, 0, 1);
		for(int i=0; i<separatrix1.size(); ++i) {
		dirEdge& edge = separatrix1[i];
		Delaunay::Vertex_handle v1 = edge.first;
		Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, edge.second);
		glVertex2f(v1->point().x(), v1->point().y());
		glVertex2f(v2->point().x(), v2->point().y());
		}*/
		glEnd();

		//
		std::cout << "separatrix length = " << separatrix0.size() << std::endl;

	}

	void DelaunayGraphics::draw_edge_stretch() {
		delaunay_->calc_edge_stretch();
		delaunay_->calc_vertex_singularity();
		delaunay_->calc_face_group_mark();
		int hmark = delaunay_->face_group_size();
		std::cout << "hmark = " << hmark << std::endl;

		//calc the most stretched edge for each face group
		std::vector<EdgeValue> isolate_edges(hmark, std::make_pair(triEdge(), 0));
		FOR_EACH_FACE(Delaunay, delaunay_, it) {
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

		//draw grouped faces
		glPolygonMode(GL_FRONT_AND_BACK ,GL_FILL);
		glBegin(GL_TRIANGLES);
		FOR_EACH_FACE(Delaunay, delaunay_, it) {
			int gmark = it->group_mark;
			if (gmark != -1 && gmark <= stretch_all_rank_)
			{

				glColor3f(gmark * 1.0/hmark, 1 - gmark*1.0/hmark, 0);
				for(int i=0; i<3; ++i) {
					Delaunay::Vertex_handle v = it->vertex(i);
					vec2 p = to_geex(v->point());
					glVertex2f(p.x, p.y);
				}
			}
		}
		glEnd();

		//collect valid stretched edges
		std::priority_queue<EdgeValue, std::vector<EdgeValue>, EdgeValue_greater> short_edges;
		std::priority_queue<EdgeValue, std::vector<EdgeValue>, EdgeValue_less> long_edges;
		FOR_EACH_FACE(Delaunay, delaunay_, it) {
			if (delaunay_->is_infinite(it)) continue;
			if (delaunay_->has_irregular(it)) {
				for (int i=0; i<3; ++i)
				{
					Delaunay::Vertex_handle v1 = it->vertex((i+1)%3);
					Delaunay::Vertex_handle v2 = it->vertex((i+2)%3);

					//
					if (v1->singularity != 0 
						|| (v1->dual_intersects_boundary && v2->dual_intersects_boundary)) {

							double stretch = it->edge_stretch[i]; 

							if (stretch > 0)
							{
								long_edges.push(std::make_pair(triEdge(it, i), stretch));
							}
							if (stretch < 0)
							{
								short_edges.push(std::make_pair(triEdge(it, i), stretch));
							}
					}
				}
			}
		}

		//Covert to sorted vector of unique edges
		///long edges
		std::vector<EdgeValue> longVec;
		if (!long_edges.empty()) {
			longVec.push_back(long_edges.top());
			long_edges.pop();
		}
		while (!long_edges.empty())
		{
			if (long_edges.top().first != longVec.back().first)
			{
				longVec.push_back(long_edges.top());
			}
			long_edges.pop();
		}
		//short edges
		std::vector<EdgeValue> shortVec;
		if(!short_edges.empty()) {
			shortVec.push_back(short_edges.top());
			short_edges.pop();
		}
		while (!short_edges.empty())
		{
			if (short_edges.top().first != shortVec.back().first)
			{
				shortVec.push_back(short_edges.top());
			}
			short_edges.pop();
		}

		//Draw Long Edges
		if (!longVec.empty() && stretch_long_rank_ >= 0)
		{
			//calc max stretch
			double max_stretch =  longVec[0].second;
			if (stretch_long_rank_ < longVec.size())
			{
				max_stretch = longVec[stretch_long_rank_].second;
				std::cout << "stretch: max[" << stretch_long_rank_ << "] " << max_stretch;
			}
			//draw
			glLineWidth(2.5) ;
			glBegin(GL_LINES);
			//draw front edges
			int front_size = (1+stretch_long_rank_ < longVec.size() ? (1+stretch_long_rank_) : longVec.size());
			for (int i=0; i<front_size; ++i)
			{
				glColor3f(0, 1, 1);
				draw_triEdge(longVec[i].first);
			}
			//draw rest edges
			for (int i=front_size; i<longVec.size(); ++i)
			{
				EdgeValue& ev = longVec[i];
				triEdge& e = ev.first;
				double stretch = ev.second;
				glColor3f(0, 0, stretch/max_stretch);
				draw_triEdge(e);
			}
			glEnd();

		}

		//Draw Short Edges
		if (!longVec.empty() && stretch_short_rank_ >= 0)
		{
			//calc min stretch
			double min_stretch =  shortVec[0].second;
			if (stretch_short_rank_ < shortVec.size())
			{
				min_stretch = shortVec[stretch_short_rank_].second;
				std::cout << " min[" << stretch_short_rank_ << "] " << min_stretch << std::endl;
			}
			//draw
			glLineWidth(2.5) ;
			glBegin(GL_LINES);
			//draw front edges
			int front_size = (1+stretch_short_rank_ < shortVec.size() ? (1+stretch_short_rank_) : shortVec.size());
			for (int i=0; i<front_size; ++i)
			{
				glColor3f(1, 0, 1);
				draw_triEdge(shortVec[i].first);
			}
			//draw rest edges
			for (int i=front_size; i<shortVec.size(); ++i)
			{
				EdgeValue& ev = shortVec[i];
				triEdge& e = ev.first;
				double stretch = ev.second;
				glColor3f(stretch/min_stretch, 0, 0);
				draw_triEdge(e);
			}
			glEnd();
		}

		//All Edges
		//merge short and long edges
		// 		std::vector<EdgeValue> allVec;
		// 		std::vector<EdgeValue>::iterator lIt = longVec.begin(), sIt = shortVec.begin();
		// 		while(true) {
		// 			if (lIt == longVec.end())
		// 			{
		// 				while (sIt != shortVec.end())
		// 				{
		// 					allVec.push_back(*sIt);
		// 					++sIt;
		// 				}
		// 				break;
		// 			}
		// 			if (sIt == shortVec.end())
		// 			{
		// 				while (lIt != longVec.end())
		// 				{
		// 					allVec.push_back(*lIt);
		// 					++lIt;
		// 				}
		// 				break;
		// 			}
		// 
		// 			if (lIt->second > abs(sIt->second))
		// 			{
		// 				allVec.push_back(*lIt);
		// 				++lIt;
		// 			}
		// 			else {
		// 				allVec.push_back(*sIt);
		// 				++sIt;
		// 			}
		// 		}
		/*if (!allVec.empty() && stretch_all_rank_ >= 0)
		{
		int front_size = (1+stretch_all_rank_ < allVec.size()) ? (1+stretch_all_rank_) : allVec.size();
		std::cout << "stretch all[" << front_size -1 << "] " << allVec[front_size-1].second << std::endl;
		glLineWidth(2.5);
		glBegin(GL_LINES);
		for(int i=0; i<front_size; ++i) {
		EdgeValue& ev = allVec[i];
		if (ev.second > 0)
		{
		glColor3f(0, 1, 0.8);
		}
		else {
		glColor3f(1, 0.5, 0);
		}
		draw_triEdge(ev.first);
		}
		for(int i=front_size; i<allVec.size(); ++i) {
		glColor3f(0, 0, 0);
		draw_triEdge(allVec[i].first);
		}
		glEnd();
		}*/

		//draw isolated stretched edges
		if (!isolate_edges.empty() || stretch_all_rank_ >= 0)
		{
			int front_size = (1+stretch_all_rank_ < isolate_edges.size()) ? (1+stretch_all_rank_) : isolate_edges.size();
			glLineWidth(2.5);
			glBegin(GL_LINES);
			for (int i=0; i<front_size; ++i)
			{
				triEdge& e = isolate_edges[i].first;
				double stretch = isolate_edges[i].second;

				if (stretch > 0)
				{
					glColor3f(0, 1, 0.8);
				} 
				else
				{
					glColor3f(1, 0.5, 0);
				}
				draw_triEdge(e);
			}
			glEnd();
		}

	}

	void DelaunayGraphics::draw_selected() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glPointSize(10); ///test
		glLineWidth(3);

		FOR_EACH_FACE(Delaunay, delaunay_, it) {
			if (delaunay_->is_infinite(it)) continue;

			//Vertex
			std::vector<Point> points;
			glColor3f(0.5, 0.0, 1.0);
			glBegin(GL_POINTS);
			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v = it->vertex(i);
				points.push_back(v->point());
				if (v->selected)
				{
					glVertex2f(v->point().x(), v->point().y());
				}
			}
			glEnd();

			//Edge
			glColor3f(0.0, 0.0, 0.0);
			glBegin(GL_LINES);
			for (int i=0; i<3; ++i)
			{
				if (it->selected_edge[i])
				{
					int j1 = (i+1) % 3;
					int j2 = (i+2) % 3;
					glVertex2f(points[j1].x(), points[j1].y());
					glVertex2f(points[j2].x(), points[j2].y());
				}

			}
			glEnd();

			//Face
			if (it->selected)
			{
				glColor3f(1.0, 0.5, 0.0);
				glBegin(GL_TRIANGLES);
				for (int i=0; i<3; ++i)
				{
					glVertex2f(points[i].x(), points[i].y());
				}
				glEnd();
			}

		}

	}

	void DelaunayGraphics::rotate_vector(vec2& v, double angle) {
		double c = cos(angle);
		double s = sin(angle);
		double new_x = c * v.x - s * v.y;
		double new_y = s * v.x + c * v.y;
		v = vec2(new_x, new_y);
	}

	void DelaunayGraphics::draw_vector(const vec2& v, const vec2& start) {
		double len = v.length();
		if (len == 0) return;
		float arror_size = 0.15;
		vec2 arror1 = v;
		rotate_vector(arror1, 0.9 * M_PI);
		arror1 *= arror_size;
		vec2 arror2 = v;
		rotate_vector(arror2, -0.9 * M_PI);
		arror2 *= arror_size;

		glBegin(GL_LINES) ;
		vec2 end = start + v;
		glVertex2f(start.x, start.y);
		glVertex2f(end.x, end.y);
		vec2 arror1_end = end + arror1;
		glVertex2f(end.x, end.y);
		glVertex2f(arror1_end.x, arror1_end.y);
		vec2 arror2_end = end + arror2;
		glVertex2f(end.x, end.y);
		glVertex2f(arror2_end.x, arror2_end.y);
		glEnd();
	}

	void DelaunayGraphics::draw_primal_triangle(Delaunay::Face_handle f) {
		glBegin(GL_TRIANGLES);
		for (int i=0; i<3; ++i)
		{
			glVertex2f(f->vertex(i)->point().x(), f->vertex(i)->point().y());
		}
		glEnd();
	}
	//triEdge
	void DelaunayGraphics::draw_triEdge(const triEdge& e) {
		Delaunay::Vertex_handle v1 = e.start();
		Delaunay::Vertex_handle v2 = e.end();
		glVertex2f(v1->point().x(), v1->point().y());
		glVertex2f(v2->point().x(), v2->point().y());
	}

	///dxy add end
}


