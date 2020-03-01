#include<stdexcept>
#include<sudodem/core/Material.hpp>
#include<sudodem/core/Scene.hpp>

const shared_ptr<Material> Material::byId(int id, Scene* w_){
	Scene* w=w_?w_:Omega::instance().getScene().get();
	assert(id>=0 && (size_t)id<w->materials.size());
	assert(w->materials[id]->id == id);
	return w->materials[id];
}

const shared_ptr<Material> Material::byLabel(const std::string& label, Scene* w_){
	Scene* w=w_?w_:Omega::instance().getScene().get();
	FOREACH(const shared_ptr<Material>& m, w->materials){
		if(m->label == label) return m;
	}
	throw std::runtime_error(("No material labeled `"+label+"'.").c_str());
}

const int Material::byLabelIndex(const std::string& label, Scene* w_){
	Scene* w=w_?w_:Omega::instance().getScene().get(); size_t iMax=w->materials.size();
	for(size_t i=0; i<iMax; i++){
		if(w->materials[i]->label==label) return i;
	}
	throw std::runtime_error(("No material labeled `"+label+"'.").c_str());
}
