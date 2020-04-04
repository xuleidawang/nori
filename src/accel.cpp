/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

bool Accel::rayIntersect(OcTreeNode* node, const Ray3f &ray_, Intersection &its, bool shadowRay) const{
    if(!node) return false;
    bool foundIntersection = false;
    uint32_t f = -1;

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)
    //Ray doesn't hit the bounding box, reject
    if(!node->bbox.rayIntersect(ray))
        return false;

    //Interior node
    if(node->triangleIndices.size() < 1)
    {
        std::vector<std::pair<float, int> > boundsHitDistance;
        for(int i=0; i< 8; i++)
        {
            Ray3f ray(ray_);
            float nearT, farT;
            // Create child node bound index and  hit distance pair
            if( node->children[i] && node->children[i]->bbox.rayIntersect(ray, nearT, farT))
            {
                boundsHitDistance.push_back(std::make_pair(nearT, i));
            }
            
        }
        //sort to traverse the nodes from near to far
        sort(boundsHitDistance.begin(), boundsHitDistance.end());

        for(auto p: boundsHitDistance)
        {
            Ray3f ray(ray_); 
            Intersection its_;
            if (rayIntersect(node->children[p.second], ray, its_, shadowRay))
            {
                its = its_;
                return true;
            } 
        }
        return false;
    }
    //It is a leaf node
    else
    {
        const Mesh* m = this->getMesh();
        Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)
        for(uint32_t i=0;i < node->triangleIndices.size(); i++)
        {    
            float u, v, t;
            if(m->rayIntersect(node->triangleIndices[i], ray, u, v, t))
            {
                if(shadowRay) return true;

                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                f = node->triangleIndices[i];
                foundIntersection = true;
            }
        }

        if (foundIntersection)
        {
            Vector3f bary;
            bary << 1-its.uv.sum(), its.uv;

            /* References to all relevant mesh buffers */
            const Mesh *mesh   = its.mesh;
            const MatrixXf &V  = mesh->getVertexPositions();
            const MatrixXf &N  = mesh->getVertexNormals();
            const MatrixXf &UV = mesh->getVertexTexCoords();
            const MatrixXu &F  = mesh->getIndices();
            /* Vertex indices of the triangle */
            uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);
            Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);
            /* Compute the intersection positon accurately
            using barycentric coordinates */
            its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;
            /* Compute proper texture coordinates if provided by the mesh */
            if (UV.size() > 0)
                its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

            /* Compute the geometry frame */
            its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

            if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

                its.shFrame = Frame(
                    (bary.x() * N.col(idx0) +
                     bary.y() * N.col(idx1) +
                    bary.z() * N.col(idx2)).normalized());
            } 
            else 
            {
                its.shFrame = its.geoFrame;
            }
        }
        
    }
    return foundIntersection;
}

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    /* Nothing to do here for now */
    std::vector<uint32_t> triangleIndices;
    for(uint32_t i=0; i < m_mesh->getTriangleCount(); i++)
    {
        triangleIndices.push_back(i);
    }
    this->root = recursiveBuild(m_mesh->getBoundingBox(), triangleIndices, 0);
    std::cout<<"OcTree Build Complete";
}
/*
Define order of sub-node bounding box.
bottom layer
| 2 | 3 |
  -   -  
| 0 | 1 |   

top layer
| 6 | 7 |
  -   -  
| 4 | 5 |   
*/

OcTreeNode* Accel::recursiveBuild(BoundingBox3f bbox, std::vector<uint32_t>& triangleIndices, int depth){
    if (triangleIndices.size() < 1)
        return nullptr;

    if (triangleIndices.size() < 10){
        OcTreeNode* node = new OcTreeNode(triangleIndices);
        return node;
    }

    if( depth > 20 ){
        OcTreeNode *node = new OcTreeNode(triangleIndices);
        return node;
    }

    BoundingBox3f box[8];
    
    //Calculate the 8 sub-nodes' bounding box
    for (int i=0 ; i < 8; i++)
    {
        box[i] = BoundingBox3f(0.5*(bbox.getCorner(0) + bbox.getCorner(i)),
            0.5*(bbox.getCorner(i) + bbox.getCorner(7)));
    }

    std::vector < std::vector<uint32_t> > triangle_list(8, std::vector<uint32_t>());

    for (uint32_t idx = 0; idx< triangleIndices.size(); idx++) {
        BoundingBox3f triangleBounds = this->getMesh()->getBoundingBox(idx);

        for (int i = 0; i < 8; ++i) {
            //if triangle bound overlap with sub-node bound
            if (box[i].overlaps(triangleBounds) || i==7)
            {
                //add triangle to that list.
                triangle_list[i].push_back(idx);
            }
        }
    }

    OcTreeNode *node = new OcTreeNode();

    for (int i = 0; i < 8; ++i)
    {
        node->children[i] = recursiveBuild(box[i], triangle_list[i], depth +1);
    }
    return node;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    return rayIntersect(root, ray_, its, shadowRay);    
}

NORI_NAMESPACE_END

