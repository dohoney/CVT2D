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
 *  If you modify this software, you should include notice giving the
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

#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <glut_viewer/tweak_bar.h>
#include "cvt.h"

void Lloyd() ;

namespace Geex {

    class CVTApp : public GeexApp {
    public:
        CVTApp(int argc, char** argv) : GeexApp(argc, argv) { 
            hdr_ = false ;
            boundary_filename_ = get_file_arg("line") ;
	    affinity_filename_ = get_file_arg("txt") ;
            if(boundary_filename_.length() > 0) {
                if(!Geex::FileSystem::is_file(boundary_filename_)) {
                    boundary_filename_ = Geex::FileSystem::get_project_root() + 
                        "/gx_cvt2d/" + boundary_filename_ ;
                }
            }
            nb_points_ = 64 ;
            get_arg("nb_pts", nb_points_) ;
            nb_iter_ = 10 ;
            get_arg("nb_iter", nb_iter_) ;            
            non_convex_ = GL_FALSE ;
            get_arg("non_convex", non_convex_) ;
            edit_ = GL_FALSE ;
			///dxy add
			select_ = GL_FALSE ;
			///
	    if(affinity_filename_.length() > 0) {
		nb_points_ = 256 ;
	    }
        }

        CVT* cvt() { return static_cast<CVT*>(scene()) ; }

        GLboolean& edit() { return edit_ ; }
		///dxy
		GLboolean& select() { return select_ ; }

        virtual void init_scene() {
            scene_ = new CVT ;
            std::cerr << "Non convex = " 
                      << (non_convex_ ? "true" : "false") 
                      << "(use +non_convex to set)" << std::endl ;
            cvt()->set_non_convex_mode(non_convex_) ;
            if(boundary_filename_.length() > 0) {
				///dxy change
                //cvt()->load_boundary(boundary_filename_) ;
				cvt()->load_boundary(boundary_filename_, nb_points_);
				///
            }
	    if(affinity_filename_.length() > 0) {
		cvt()->load_affinity(affinity_filename_) ;
	    }
            cvt()->insert_random_vertices(nb_points_) ;
			///dxy comment: done in end_insert()
			//cvt()->clear_select();
			///
        }

        void Lloyd() {
            cvt()->lloyd(nb_iter_) ;
        }

        void Lloyd(int nb_iter) {
            cvt()->lloyd(nb_iter) ;
        }

		///dxy add
		void Update_Energy_Grad() {
			std::cout << "Energy = " << cvt()->lloyd_energy() << std::endl ;
		}

		void Update_Total_Area() {
			cvt()->delaunay()->update_area();
		}

		int& nb_iter() { return nb_iter_ ; }
		///

        void NewtonLloyd() {
			cvt()->lloyd_energy(); //trick: otherwise, update_area() will not work
			cvt()->delaunay()->update_area();
			if (cvt()->use_topo_optimization())
			{
				cvt()->newton_lloyd(1);
				int it = 0;
				while (it < nb_iter_)
				{
					cvt()->stretch_topo_optimize();
					cvt()->newton_lloyd(5);
					it += 5;
				}
			}
			else {
				cvt()->newton_lloyd(nb_iter_) ;
			}
        }

        void reset() {
            cvt()->clear() ;
            cvt()->insert_random_vertices(nb_points_) ;
		}

        void grid() {
            cvt()->clear() ;
            cvt()->insert_grid() ;
        }

        virtual void init_gui() {
            GeexApp::init_gui() ;

            // New-style GUI =====================================================================================

            TwBar* graphics_bar = TwNewBar("Graphics") ;
            TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_FLOAT, &cvt()->vertices_size(), "min=0 max=1 step=0.01") ;
            TwAddVarRW(graphics_bar, "Centers", TW_TYPE_FLOAT, &cvt()->centers_size(), "min=0 max=1 step=0.01") ;
            TwAddVarRW(graphics_bar, "Primal", TW_TYPE_BOOL8, &cvt()->show_primal_mesh(), "") ;
            TwAddVarRW(graphics_bar, "Dual", TW_TYPE_BOOL8, &cvt()->show_dual_mesh(), "") ;
			TwAddVarRW(graphics_bar, "Non-Hex", TW_TYPE_BOOL8, &cvt()->show_cells(), "");  ///dxy restore
            TwAddVarRW(graphics_bar, "Energy", TW_TYPE_BOOL8, &cvt()->show_energy(), "") ;
			///dxy add
			TwAddVarRW(graphics_bar, "Area", TW_TYPE_BOOL8, &cvt()->show_relative_area(), "");  
			TwAddVarRW(graphics_bar, "Tri Area", TW_TYPE_BOOL8, &cvt()->show_triangle_area(), ""); 
			///
			///dxy comment
            /*TwAddVarRW(graphics_bar, "Bkgnd. field", TW_TYPE_BOOL8, &cvt()->show_field(), "") ;
            TwAddVarRW(graphics_bar, "Quads", TW_TYPE_FLOAT, &cvt()->quad_ratio(), "min=0.0 max=1.5 step=0.01") ;
            TwAddVarRW(graphics_bar, "show Lp", TW_TYPE_BOOL8, &cvt()->show_lp(), "") ;            
            TwAddVarRW(graphics_bar, "Lp view", TW_TYPE_INT32, &cvt()->lp_shader(), "min=0 max=10 step=1") ;
            TwAddVarRW(graphics_bar, "Lp scale", TW_TYPE_FLOAT, &cvt()->new_uniform("scale", 2.0), "min=0.001 max=10 step=0.1") ;
            TwAddVarRW(graphics_bar, "show affinity", TW_TYPE_BOOL8, &cvt()->show_affinity(), "") ;         */
			///dxy add
			TwAddSeparator(graphics_bar, "Grad", "");
			TwAddVarRW(graphics_bar, "Lloyd Grad Mag", TW_TYPE_INT32, &cvt()->lloyd_grad_magnification(), "min=1 max=100 step=1");
			TwAddVarRW(graphics_bar, "Direct Grad Mag", TW_TYPE_INT32, &cvt()->direct_grad_magnification(), "min=1 max=100 step=1");
			TwAddVarRW(graphics_bar, "Total Grad Mag", TW_TYPE_INT32, &cvt()->total_grad_magnification(), "min=1 max=100 step=1");
			TwAddVarRW(graphics_bar, "Lloyd Gradient", TW_TYPE_BOOL8, &cvt()->show_lloyd_grad(), "") ;
			TwAddVarRW(graphics_bar, "Direction Gradient", TW_TYPE_BOOL8, &cvt()->show_direction_grad(), "") ;
			TwAddVarRW(graphics_bar, "Total Gradient", TW_TYPE_BOOL8, &cvt()->show_total_grad(), "");
			TwAddSeparator(graphics_bar, "Edge", "");
			TwAddVarRW(graphics_bar, "Separatrix", TW_TYPE_BOOL8, &cvt()->show_separatrix(), "");
			TwAddVarRW(graphics_bar, "Edge Stretch", TW_TYPE_BOOL8, &cvt()->show_edge_stretch(), "");
			TwAddVarRW(graphics_bar, "- long  rank", TW_TYPE_INT32, &cvt()->stretch_long_rank(), "min=-1 max=100 step=1");
			TwAddVarRW(graphics_bar, "- short rank", TW_TYPE_INT32, &cvt()->stretch_short_rank(), "min=-1 max=100 step=1");
			TwAddVarRW(graphics_bar, "- all   rank", TW_TYPE_INT32, &cvt()->stretch_all_rank(), "min=-1 max=100 step=1");

			TwAddVarRW(graphics_bar, "show Constrained Edge", TW_TYPE_BOOL8, &cvt()->show_constrained_edge(), "");
			TwAddVarRW(graphics_bar, "show Direction Field", TW_TYPE_BOOL8, &cvt()->show_direction_field(), "");
			///

			/* TwBar* numerics_bar = TwNewBar("Numerics") ;
			TwDefine(" Numerics position='40 400' size='200 200'") ;            
			TwAddVarRW(numerics_bar, "Lp order", TW_TYPE_FLOAT, &cvt()->Lp(), "min=0 max=9 step=1") ;
			TwEnumVal aniso_def[] = { {CONSTANT, "const"}, {R_INV, "1/R"}, {BORDER, "Border"} } ;
			TwType tw_aniso = TwDefineEnum("AnisoType", aniso_def, 3) ;
			TwAddVarRW(numerics_bar, "Aniso", tw_aniso, &cvt()->aniso_mode(), "") ;
			TwType tw_center_mode = TwDefineEnum("CenterMode", "centroid, quasi-incenter") ;
			TwAddVarRW(numerics_bar, "CenterMode", tw_center_mode, &cvt()->center_mode(), "") ;
			TwAddVarRW(numerics_bar, "X scale", TW_TYPE_FLOAT, &cvt()->Xscale(), "min=-2 max=2 step=1") ;            
			TwAddVarRW(numerics_bar, "Y scale", TW_TYPE_FLOAT, &cvt()->Yscale(), "min=-2 max=2 step=1") ;            
			TwEnumVal cvt_def[] = { {LLOYD, "Lloyd"}, {NEWTON, "Newton"}, {THETA, "Theta"} } ;
			TwType tw_cvt = TwDefineEnum("CVTType", cvt_def, 3) ;
			TwAddVarRW(numerics_bar, "Mode", tw_cvt, &cvt()->mode(), "") ;
			TwAddVarRW(numerics_bar, "Bndry.", TW_TYPE_BOOL8, &cvt()->insert_boundary(), "") ;            
			TwAddVarRW(numerics_bar, "Affinity W", TW_TYPE_FLOAT, &cvt()->affinity_weight(), "min=0.0 max=10000.0 step=0.0001") ;*/
		

			///dxy add: direction bar
			TwBar* direction_bar = TwNewBar("Direction Field") ;
			TwAddVarRW(direction_bar, "nb iter", TW_TYPE_INT32, &nb_iter(), "min=0 max=1000 step=1");
			TwAddVarRW(direction_bar, "Lloyd Energy", TW_TYPE_BOOL8, &cvt()->use_lloyd_energy(), "");
			TwAddVarRW(direction_bar, "Direction Field", TW_TYPE_BOOL8, &cvt()->direction_field(), "");
			TwAddVarRW(direction_bar, "True Gradient", TW_TYPE_BOOL8, &cvt()->use_true_grad(), "");
			TwType tw_direction_field_mode = TwDefineEnum("FieldMode", "constant, circle, sin");
			TwAddVarRW(direction_bar, "Field Type", tw_direction_field_mode, &cvt()->direction_field_mode(), "");
			TwAddVarRW(direction_bar, "Direction Angle", TW_TYPE_FLOAT, &cvt()->direction_angle(), "min=-1.0 max=1.0 step=0.1");
			TwAddVarRW(direction_bar, "Direction Factor", TW_TYPE_FLOAT, &cvt()->direction_factor(), "min=0.0 max=10.0 step=0.0001");
			TwType tw_direction_energy_mode = TwDefineEnum("EnergyMode", "angle, angle-dual, angle-primal");
			TwAddVarRW(direction_bar, "EnergyMode", tw_direction_energy_mode, &cvt()->direction_energy_mode(), "");
			TwAddVarRW(direction_bar, "Direction Gradient", TW_TYPE_BOOL8, &cvt()->direction_gradient(), "");
			TwAddVarRW(direction_bar, "Constrained Edge", TW_TYPE_BOOL8, &cvt()->direction_constrained_edge(), "");
			TwAddVarRW(direction_bar, "Boundary Edge", TW_TYPE_BOOL8, &cvt()->direction_boundary_edge(), "");
			///
			TwAddSeparator(direction_bar, "Quad", "");
			TwAddVarRW(direction_bar, "Quad Field", TW_TYPE_BOOL8, &cvt()->use_quad_field(), "");
			TwAddVarRW(direction_bar, "Quad Factor", TW_TYPE_FLOAT, &cvt()->quad_factor(), "min=0.0 max=1.0 step=0.1");
			///
			TwAddSeparator(direction_bar, "Topo", "");
			TwAddVarRW(direction_bar, "Topo Optimization", TW_TYPE_BOOL8, &cvt()->use_topo_optimization(), "");

			///dxy add: edit bar
			TwBar* edit_bar = TwNewBar("Edit Topology");
			TwAddVarRW(edit_bar, "Edit", TW_TYPE_BOOL8, &edit_, "") ;  //move to here
			TwAddVarRW(edit_bar, "select", TW_TYPE_BOOL8, &select_, ""); //like edit_
			TwType tw_select_mode = TwDefineEnum("SelectMode", "vertex, edge, face, separatrix");
			TwAddVarRW(edit_bar, "SelectMode", tw_select_mode, &cvt()->delaunay()->select_mode(), "");
			TwAddVarRW(edit_bar, "show selected", TW_TYPE_BOOL8, &cvt()->show_selected(), "");
			/*TwAddVarRW(edit_bar, "use color", TW_TYPE_BOOL8, &cvt()->use_select_color(), "");
			TwAddVarRW(edit_bar, "R", TW_TYPE_INT32, &cvt()->select_color_R(), "min=0 max=255 step=1");
			TwAddVarRW(edit_bar, "G", TW_TYPE_INT32, &cvt()->select_color_G(), "min=0 max=255 step=1");
			TwAddVarRW(edit_bar, "B", TW_TYPE_INT32, &cvt()->select_color_B(), "min=0 max=255 step=1");*/
			TwAddSeparator(edit_bar, "Short_Long", "");
			TwAddVarRW(edit_bar, "Long Edge", TW_TYPE_BOOL8, &cvt()->show_long_edge(), ""); 
			TwAddVarRW(edit_bar, "Short Edge", TW_TYPE_BOOL8, &cvt()->show_short_edge(), ""); 

            viewer_properties_->add_separator("Graphics") ;
//            viewer_properties_->add_toggle("Domain", cvt()->show_domain()) ;
            viewer_properties_->add_slider("Vertices", cvt()->vertices_size()) ;
            viewer_properties_->add_slider("Centers", cvt()->centers_size()) ;
            viewer_properties_->add_toggle("Primal mesh", cvt()->show_primal_mesh()) ;
            viewer_properties_->add_toggle("Dual mesh", cvt()->show_dual_mesh()) ;
            viewer_properties_->add_toggle("Non-hex", cvt()->show_cells()) ;  ///dxy restore
            viewer_properties_->add_toggle("Energy", cvt()->show_energy()) ;
			viewer_properties_->add_toggle("Area", cvt()->show_relative_area()); ///dxy add
			viewer_properties_->add_toggle("Tri Area", cvt()->show_triangle_area()); //dxy add
            viewer_properties_->add_toggle("Bkgnd. field", cvt()->show_field()) ;
            viewer_properties_->add_slider("Quads", cvt()->quad_ratio(), 0.0, 1.5) ;
			///dxy add
			/*viewer_properties_->add_toggle("Lloyd Gradient", cvt()->show_lloyd_grad()) ;
			viewer_properties_->add_toggle("Direction Gradient", cvt()->show_direction_grad());
			viewer_properties_->add_toggle("show Constrained Edge", cvt()->show_constrained_edge());
			viewer_properties_->add_toggle("show Direction Field", cvt()->show_direction_field());*/
			///

            viewer_properties_->add_separator("Optimizer") ;
            viewer_properties_->add_slider("Lp order", cvt()->Lp(), 0.0, double(DelaunayCVT::MAX_P))->set_integer(GL_TRUE) ;
            viewer_properties_->add_enum("Aniso", cvt()->aniso_mode(), GlutViewerGUI::LabelList() | AnisoModeNames) ;            
            viewer_properties_->add_slider("X scale", cvt()->Xscale(), -2.0, 2.0)->set_integer(GL_TRUE) ;
            viewer_properties_->add_slider("Y scale", cvt()->Yscale(), -2.0, 2.0)->set_integer(GL_TRUE) ;
            viewer_properties_->add_enum("Mode", cvt()->mode(), GlutViewerGUI::LabelList() | CVTModeNames) ;
            viewer_properties_->add_toggle("Bndry.", cvt()->insert_boundary()) ;
            viewer_properties_->add_toggle("Edit", edit_) ;
            viewer_properties_->add_toggle("Affinity", cvt()->show_affinity()) ;

			///dxy add
			/*viewer_properties_->add_toggle("Lloyd Energy", cvt()->use_lloyd_energy());
			viewer_properties_->add_toggle("Direction Field", cvt()->direction_field());
			viewer_properties_->add_slider("Direction Angle", cvt()->direction_angle());
			viewer_properties_->add_slider("Direction Factor", cvt()->direction_factor());
			viewer_properties_->add_toggle("Constrained Edge", cvt()->direction_constrained_edge());
			viewer_properties_->add_toggle("Boundary Edge", cvt()->direction_boundary_edge());
			viewer_properties_->add_toggle("Direction Gradient", cvt()->direction_gradient());*/
			///

			///dxy add: quad
			/*viewer_properties_->add_toggle("Quad Field", cvt()->use_quad_field());
			viewer_properties_->add_slider("Quad Factor", cvt()->quad_factor());*/
			///

            toggle_skybox_CB() ;

            glut_viewer_add_toggle('B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
//            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
            viewer_properties_->hide() ;
            glut_viewer_add_toggle('T', glut_viewer_is_enabled_ptr(GLUT_VIEWER_TWEAKBARS), "Tweak bars") ;

            glut_viewer_disable(GLUT_VIEWER_BACKGROUND) ;
        }


        virtual GLboolean mouse(float x, float y, int button, enum GlutViewerEvent event) {
            static int timestamp = 0 ;
            static int last_timestamp = 0 ;
            timestamp++ ;
            static int mode = 0 ;


            GLdouble p[3] ;
            GLdouble v[3] ;
            if(GeexApp::mouse(x, y, button, event)) { return GL_TRUE ; }
            if(edit_) {
                
                bool change = false ;
                if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; }
                
                if(event == GLUT_VIEWER_UP) { return GL_FALSE ; }

                glut_viewer_get_picked_ray(p,v) ;
                vec2 pt(p[0], p[1]) ;

                if(mode == 0 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
                        cvt()->begin_insert() ;
                        cvt()->insert(pt) ;
                        cvt()->end_insert() ;
                    }
                    last_timestamp = timestamp ;
                }

                if(mode == 1 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {


//                        cvt()->nearest(pt)->theta += M_PI / 4.0 ;


                        cvt()->begin_insert() ;
                        cvt()->remove(pt) ;
                        cvt()->end_insert(false) ;
                        cvt()->begin_insert() ;
                        cvt()->insert(pt) ; // ->locked = true ;
                        cvt()->end_insert() ;


                    }
                    last_timestamp = timestamp ;
                }

                if(mode == 2 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
                        cvt()->begin_insert() ;
                        cvt()->remove(pt) ;
                        cvt()->end_insert() ;
                    }
                    last_timestamp = timestamp ;
                }


                return GL_TRUE ;
            }
			///dxy add: select
			if(!edit_ && select_) {

				bool change = false ;
				if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; }

				if(event == GLUT_VIEWER_UP) { return GL_FALSE ; }

				glut_viewer_get_picked_ray(p,v) ;
				vec2 pt(p[0], p[1]) ;

				if(mode == 0 && (change || timestamp - last_timestamp > 3)) { //left click - select/unselect
					if(cvt()->in_boundary(pt)) {
						cvt()->update_select(pt);
					}
					last_timestamp = timestamp ;
				}



				return GL_TRUE ;
			}
			///dxy add end
            return GL_FALSE ;
        }

    private:
        std::string boundary_filename_ ;
        int nb_points_ ;
        int nb_iter_ ;
        GLboolean non_convex_ ;
        GLboolean edit_ ;
		///dxy add
		GLboolean select_;
		///
	std::string affinity_filename_ ;
    } ;
}

Geex::CVTApp* cvt_app() { return static_cast<Geex::CVTApp*>(Geex::GeexApp::instance()) ; }

void Lloyd() {
    cvt_app()->Lloyd() ;
}

void Lloyd1() {
    cvt_app()->Lloyd(1) ;
}

///dxy add
void Update_Energy_Grad() {
	cvt_app()->Update_Energy_Grad();
}

void Update_Total_Area() {
	cvt_app()->Update_Total_Area();
}

void edge_split() {
	if (cvt_app()->cvt()->show_long_edge())
	{
		//cvt_app()->cvt()->split_long_edge();
		cvt_app()->cvt()->split_isolate_long_edge();
		return;
	}
	cvt_app()->cvt()->split_selected_edge();
	
}

void edge_collapse() {
	if (cvt_app()->cvt()->show_short_edge())
	{
		//cvt_app()->cvt()->collapse_short_edge();
		cvt_app()->cvt()->collapse_isolate_short_edge();
		return;
	}
	cvt_app()->cvt()->collapse_selected_edge();
	
}

void edge_flip() {
	cvt_app()->cvt()->flip_selected_edge();
}
///

void NewtonLloyd() {
    cvt_app()->NewtonLloyd() ;
}

void reset() {
    cvt_app()->reset() ;
}

void inc_Lp() {
    cvt_app()->cvt()->Lp() += 1.0 ;
    cvt_app()->cvt()->Lp() = Geex::gx_min(cvt_app()->cvt()->Lp(), float(Geex::DelaunayCVT::MAX_P)) ;
    cvt_app()->cvt()->lp_shader() = int(cvt_app()->cvt()->Lp()) ;
}

void dec_Lp() {
    cvt_app()->cvt()->Lp() -= 1.0 ;
    cvt_app()->cvt()->Lp() = Geex::gx_max(cvt_app()->cvt()->Lp(), 0.0f) ;
    cvt_app()->cvt()->lp_shader() = int(cvt_app()->cvt()->Lp()) ;
}

void save() {
    cvt_app()->cvt()->save("out.pts") ;
}

void grid() {
    cvt_app()->grid() ;
}

void toggle_vertices() {
    if(cvt_app()->cvt()->vertices_size() == 0.0) {
	cvt_app()->cvt()->vertices_size() = 0.5 ;
    } else {
	cvt_app()->cvt()->vertices_size() = 0.0 ;
    }
}

void toggle_centers() {
    if(cvt_app()->cvt()->centers_size() == 0.0) {
	cvt_app()->cvt()->centers_size() = 0.8 ;
    } else {
	cvt_app()->cvt()->centers_size() = 0.0 ;
    }
}

void toggle_Newton_mode() {
    if(cvt_app()->cvt()->mode() == Geex::LLOYD) {
	cvt_app()->cvt()->mode() = Geex::NEWTON ;
    } else {
	cvt_app()->cvt()->mode() = Geex::LLOYD ;
    }
}

void toggle_Lp_display() {
    if(cvt_app()->cvt()->show_lp()) {
	cvt_app()->cvt()->show_lp() = false ;
    } else {
	cvt_app()->cvt()->show_lp() = true ;
    }
}

void toggle_primal_display() {
    if(cvt_app()->cvt()->show_primal_mesh()) {
	cvt_app()->cvt()->show_primal_mesh() = false ;
    } else {
	cvt_app()->cvt()->show_primal_mesh() = true ;
    }
}

int main(int argc, char** argv) {
    Geex::initialize() ;
    Geex::CVTApp app(argc, argv) ;
    glut_viewer_add_key_func('k', Lloyd, "Lloyd iterations") ;
    glut_viewer_add_key_func('K', Lloyd1, "Lloyd one iteration") ;
	///dxy add
	glut_viewer_add_key_func('n', Update_Energy_Grad, "Update energy and grad");
	glut_viewer_add_key_func('A', Update_Total_Area, "Update total area");
	///
    glut_viewer_add_key_func('m', NewtonLloyd, "Newton-Lloyd iterations") ;
    glut_viewer_add_key_func('Z', reset, "reset") ;
    glut_viewer_add_key_func('G', grid, "grid") ;
    glut_viewer_add_key_func('S',  save, "save") ;
    glut_viewer_add_toggle('e', &(cvt_app()->edit()), "edit") ;
    glut_viewer_add_key_func('<', dec_Lp, "decrement Lp") ;
    glut_viewer_add_key_func('>', inc_Lp, "increment Lp") ;
    glut_viewer_add_key_func('v', toggle_vertices, "toggle vertices") ;
    glut_viewer_add_key_func('g', toggle_centers, "toggle centers"); 
    glut_viewer_add_key_func('N', toggle_Newton_mode, "toggle Newton mode"); 
    glut_viewer_add_key_func('L', toggle_Lp_display, "toggle Lp display"); 
    glut_viewer_add_key_func('P', toggle_primal_display, "toggle primal display");
	///dxy add
	glut_viewer_add_key_func('s', edge_split, "edge split");
	glut_viewer_add_key_func('c', edge_collapse, "edge collapse");
	glut_viewer_add_key_func('f', edge_flip, "edge flip");
	///

    glut_viewer_disable(GLUT_VIEWER_3D) ;
    app.main_loop() ;
    Geex::terminate() ;
}
