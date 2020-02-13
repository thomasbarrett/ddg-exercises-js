"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                for (let i = 0; i < mesh.vertices.length; i++) {
                        mesh.vertices[i].index = i;
                }
                for (let i = 0; i < mesh.edges.length; i++) {
                        mesh.edges[i].index = i;
                }
                for (let i = 0; i < mesh.faces.length; i++) {
                        mesh.faces[i].index = i;
                }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                let T = new Triplet(mesh.edges.length, mesh.vertices.length);
                for (let e of mesh.edges) {
                        let v1 = e.halfedge.vertex;
                        let v2 = e.halfedge.twin.vertex;

                        T.addEntry(1, e.index, v1.index);
                        T.addEntry(1, e.index, v2.index);
                }
                return SparseMatrix.fromTriplet(T);
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let T = new Triplet(mesh.faces.length, mesh.edges.length);
                for (let f of mesh.faces) {
                        let e1 = f.halfedge.edge;
                        let e2 = f.halfedge.next.edge;
                        let e3 = f.halfedge.next.next.edge;

                        T.addEntry(1, f.index, e1.index);
                        T.addEntry(1, f.index, e2.index);
                        T.addEntry(1, f.index, e3.index);
                }
                return SparseMatrix.fromTriplet(T);
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let vec = DenseMatrix.zeros(this.mesh.vertices.length, 1);
                for (let v of subset.vertices) {
                        vec.set(1, v, 0);
                }
                return vec;
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let vec = DenseMatrix.zeros(this.mesh.edges.length, 1);
                for (let e of subset.edges) {
                        vec.set(1, e, 0);
                }
                return vec;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let vec = DenseMatrix.zeros(this.mesh.faces.length, 1);
                for (let f of subset.faces) {
                        vec.set(1, f, 0);
                }
                return vec;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // TODO

                return subset; // placeholder
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // TODO

                return subset; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // TODO

                return subset; // placeholder
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                // TODO
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                // TODO
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                let boundarySubset = new MeshSubset();
                let faceEdges = this.A1.transpose().timesDense(this.buildFaceVector(subset));
                for (let i = 0; i < this.mesh.edges.length; i++) {
                        if (faceEdges[i] == 1) {
                                boundarySubset.addEdge(i);
                        }
                }
                return boundarySubset; // placeholder
        }
}
