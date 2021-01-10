////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "ExternalModule.h"
#include <algorithm>
#if defined(SUANPAN_WIN)
#include <Windows.h>
#elif defined(SUANPAN_UNIX)
#include <dlfcn.h>
#endif

using element_creator = void (*)(unique_ptr<Element>&, istringstream&);
using load_creator = void (*)(unique_ptr<Load>&, istringstream&);
using material_creator = void (*)(unique_ptr<Material>&, istringstream&);
using section_creator = void (*)(unique_ptr<Section>&, istringstream&);
using solver_creator = void (*)(unique_ptr<Solver>&, istringstream&);
using amplitude_creator = void (*)(unique_ptr<Amplitude>&, istringstream&);
using modifier_creator = void (*)(unique_ptr<Modifier>&, istringstream&);
using constraint_creator = void (*)(unique_ptr<Constraint>&, istringstream&);

ExternalModule::ExternalModule(string L)
	: library_name(std::move(L)) {
#ifdef SUANPAN_WIN
	auto file_name = library_name + ".dll";
	auto gnu_name = "lib" + file_name;

	// no prefix
	ext_library = LoadLibraryA(file_name.c_str());
	if(ext_library == nullptr) {
		transform(file_name.begin(), file_name.end(), file_name.begin(), suanpan::to_lower);
		ext_library = LoadLibraryA(file_name.c_str());
	}
	if(ext_library == nullptr) {
		transform(file_name.begin(), file_name.end(), file_name.begin(), suanpan::to_upper);
		ext_library = LoadLibraryA(file_name.c_str());
	}

	// linux prefix style
	if(ext_library == nullptr) ext_library = LoadLibraryA(gnu_name.c_str());
	if(ext_library == nullptr) {
		transform(gnu_name.begin(), gnu_name.end(), gnu_name.begin(), suanpan::to_lower);
		ext_library = LoadLibraryA(gnu_name.c_str());
	}
	if(ext_library == nullptr) {
		transform(gnu_name.begin(), gnu_name.end(), gnu_name.begin(), suanpan::to_upper);
		ext_library = LoadLibraryA(gnu_name.c_str());
	}
	if(ext_library == nullptr) suanpan_error("locate_module() cannot find the library with the given name %s.\n", file_name.c_str());
#elif defined(SUANPAN_UNIX)
	auto file_name = "./lib" + library_name + ".so";
	ext_library = dlopen(file_name.c_str(), RTLD_NOW);
	if(ext_library == nullptr) {
		file_name = "./" + library_name + ".so";
		ext_library = dlopen(file_name.c_str(), RTLD_NOW);
	}
	if(ext_library == nullptr) suanpan_error("locate_module() cannot find the library with the given name %s.\n", file_name.c_str());
#endif
}

ExternalModule::~ExternalModule() {
#ifdef SUANPAN_WIN
	if(ext_library != nullptr) FreeLibrary(HINSTANCE(ext_library));
#elif defined(SUANPAN_UNIX)
    if(ext_library != nullptr) dlclose(ext_library);
#endif
}

bool ExternalModule::locate_module(string module_name) {
	if(ext_library == nullptr) return false;

	transform(module_name.begin(), module_name.end(), module_name.begin(), suanpan::to_lower);
	module_name = "new_" + module_name;

#ifdef SUANPAN_WIN
	ext_creator = reinterpret_cast<void*>(GetProcAddress(HINSTANCE(ext_library), LPCSTR(module_name.c_str())));
#elif defined(SUANPAN_UNIX)
    ext_creator = dlsym(ext_library, module_name.c_str());
#endif

	return ext_creator != nullptr;
}

void ExternalModule::new_object(unique_ptr<Element>& return_obj, istringstream& command) const { (element_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Load>& return_obj, istringstream& command) const { (load_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Material>& return_obj, istringstream& command) const { (material_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Section>& return_obj, istringstream& command) const { (section_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Solver>& return_obj, istringstream& command) const { (solver_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Amplitude>& return_obj, istringstream& command) const { (amplitude_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Modifier>& return_obj, istringstream& command) const { (modifier_creator(ext_creator))(return_obj, command); }

void ExternalModule::new_object(unique_ptr<Constraint>& return_obj, istringstream& command) const { (constraint_creator(ext_creator))(return_obj, command); }
