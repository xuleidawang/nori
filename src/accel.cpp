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
        triangleIndices[i] = i;
    }
    this->root = recursiveBuild(m_mesh->getBoundingBox(), triangleIndices);
}
/*
Define order of sub-node bounding box.
bottom layer
| 3 | 2 |
  -   -  
| 0 | 1 |   

top layer
| 7 | 6 |
  -   -  
| 4 | 5 |   
*/

OctreeNode* Accel::recursiveBuild(BoundingBox3f bbox, std::vector<uint32_t>& triangleIndices){
    if (triangleIndices.size() < 1)
        return nullptr;

    if (triangleIndices.size() < 10){
        OctreeNode* node = new OctreeNode(triangleIndices);
        return node;
    }

    BoundingBox3f box[8] = {bbox};
    
    for (int i=0 ; i < 8; i++)
    {
        box[i].min = 0.5*(bbox.getCorner(0) + bbox.getCorner(i));
        box[i].max = 0.5*(bbox.getCorner(i) + bbox.getCorner(7));
    }

    std::vector < std::vector<uint32_t> > triangle_list(8, std::vector<uint32_t>());
    auto m_F = this->getMesh()->getIndices();
    auto m_V = this->getMesh()->getVertexPositions();
    for (uint32_t idx = 0; idx< triangleIndices.size(); idx++) {
        for (int i = 0; i < 8; ++i) {
            //if triangle bound overlap with sub-node bound
            uint32_t i0 = m_F(0, idx), i1 = m_F(1, idx), i2 = m_F(2, idx);
            const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
            BoundingBox3f triangleBounds(p0);
            triangleBounds.expandBy(p1);
            triangleBounds.expandBy(p2);

            if (box[i].overlaps(triangleBounds))
            {
                //add triangle to that list.
                triangle_list[i].push_back(idx);
            }
        }
    }

    OctreeNode *node = new OctreeNode();

    for (int i = 0; i < 8; ++i)
        node->children[i] = recursiveBuild(box[i], triangle_list[i]);
    return node;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Brute force search through all triangles */
    for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
        float u, v, t;
        if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
            /* An intersection was found! Can terminate
               immediately if this is a shadow ray query */
            if (shadowRay)
                return true;
            ray.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = m_mesh;
            f = idx;
            foundIntersection = true;
        }
    }

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
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
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

