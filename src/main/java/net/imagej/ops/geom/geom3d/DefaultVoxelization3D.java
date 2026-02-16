/*
 * #%L
 * ImageJ2 software for multidimensional image processing and analysis.
 * %%
 * Copyright (C) 2014 - 2024 ImageJ2 developers.
 * %%
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 * #L%
 */
package net.imagej.ops.geom.geom3d;

import net.imagej.mesh.Mesh;
import net.imagej.mesh.Meshes;
import net.imagej.mesh.Triangle;
import net.imagej.ops.OpService;
import net.imagej.ops.Ops;
import net.imagej.ops.special.function.AbstractUnaryFunctionOp;
import net.imglib2.*;
import net.imglib2.img.Img;
import net.imglib2.iterator.LocalizingIntervalIterator;
import net.imglib2.type.logic.BitType;

import net.imglib2.util.Intervals;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.scijava.ItemIO;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * <p>
 * This is a voxelizer that produces a binary image with values filled in along
 * the surface of the mesh.
 * </p>
 * 
 * @author Andrew McCall (University at Buffalo)
 */
@Plugin(type = Ops.Geometric.Voxelization.class)
public class DefaultVoxelization3D extends AbstractUnaryFunctionOp<Mesh, RandomAccessibleInterval<BitType>>
		implements Ops.Geometric.Voxelization {

	@Parameter(type = ItemIO.INPUT, required = false)
	private Interval dimensions;

	@Parameter(type = ItemIO.INPUT, required = false)
	private boolean scaleMeshToDimesions = false;

	@Parameter
	private OpService ops;

	private final double wallThickness = 1.0; //This could be made into a parameter if needed
	private final long[] offset = new long[] {0,0,0};
	private double scale = 1.0;

	@Override
	public RandomAccessibleInterval<BitType> calculate(Mesh input) {

		if(dimensions == null) {
			float[] bounds = Meshes.boundingBox(input);
			long[] outputInterval = new long[3];
			for (int i = 0; i < 3; i++) {
				outputInterval[i] = (long)Math.ceil(bounds[i+3]-bounds[i]);
			}
			dimensions = new FinalInterval(outputInterval);
			setScale(input, bounds);
			scaleMeshToDimesions = false;
		}

		if(scaleMeshToDimesions)
			setScale(input);

		Img<BitType> outImg = ops.create().img(dimensions, new BitType());

		RandomAccess<BitType> ra = outImg.randomAccess();

		input.triangles().forEach((Triangle t) -> {
			Vector3D[] scaledT = scaleTriangleToOutput(t);
			Interval triangleBox = boundingBox(scaledT);
			LocalizingIntervalIterator it = new LocalizingIntervalIterator(triangleBox);
			while (it.hasNext()) {
				it.fwd();
				if(Intervals.contains(dimensions, it.positionAsPoint())) {
					if (pointToTriangleDist(new Vector3D(it.getDoublePosition(0), it.getDoublePosition(1), it.getDoublePosition(2)), scaledT) < wallThickness/2) {
						synchronized (ra) {
							ra.setPositionAndGet(it.positionAsPoint()).set(true);
						}
					}
				}
			}
		});

		return outImg;
	}

	private Vector3D[] scaleTriangleToOutput(Triangle t){
		Vector3D[] o = new Vector3D[3];

		o[0] = new Vector3D((t.v0x()-offset[0])*scale, (t.v0y()-offset[1])*scale, (t.v0z()-offset[2])*scale);
		o[1] = new Vector3D((t.v1x()-offset[0])*scale, (t.v1y()-offset[1])*scale, (t.v1z()-offset[2])*scale);
		o[2] = new Vector3D((t.v2x()-offset[0])*scale, (t.v2y()-offset[1])*scale, (t.v2z()-offset[2])*scale);
		return o;
	}

	private void setScale(Mesh input){
		setScale(input, Meshes.boundingBox(input));
	}

	private void setScale(Mesh input, float[] bounds) {
		double[] axisScaling = new double[3];
		for (int i = 0; i < 3; i++) {
			offset[i] = (long) Math.floor(bounds[i]) - dimensions.min(i);
			axisScaling[i] = (dimensions.max(i) - dimensions.min(i)) / ((bounds[i + 3]+2) - bounds[i]);
		}
		scale = Math.min(axisScaling[0], Math.min(axisScaling[1], axisScaling[2]));
	}


	private Interval boundingBox(Vector3D[] t){
		long [] min = new long[3];
		long [] max = new long[3];

		min[0] = (long) Math.floor(Math.min(t[0].getX(), Math.min(t[1].getX(), t[2].getX())));
		min[1] = (long) Math.floor(Math.min(t[0].getY(), Math.min(t[1].getY(), t[2].getY())));
		min[2] = (long) Math.floor(Math.min(t[0].getZ(), Math.min(t[1].getZ(), t[2].getZ())));

		max[0] = (long) Math.ceil(Math.max(t[0].getX(), Math.max(t[1].getX(), t[2].getX())));
		max[1] = (long) Math.ceil(Math.max(t[0].getY(), Math.max(t[1].getY(), t[2].getY())));
		max[2] = (long) Math.ceil(Math.max(t[0].getZ(), Math.max(t[1].getZ(), t[2].getZ())));
		return new FinalInterval(min, max);

	}

	private double pointToTriangleDist(Vector3D p, Vector3D[] t){
		Vector3D tPoint = nearestPointInTriangle3D(p, t);
		return p.distance(tPoint);
	}



	private Vector3D nearestPointInTriangle3D(Vector3D p, Vector3D[] t) {

		Vector3D ab = t[1].subtract(t[0]);
		Vector3D ac = t[2].subtract(t[0]);


		//region Obtain projection (p) of origP onto plane of triangle
		// Find the normal to the plane: n = (b - a) x (c - a)
		Vector3D n = ab.crossProduct(ac);

		// Normalize normal vector
		try{
			n = n.normalize();
		}
		catch(Exception e){
			return  new Vector3D(-100,-100,-100);  // Triangle is degenerate
		}

		//    Project point origP onto the plane spanned by a->b and a->c.
		double dist = p.dotProduct(n) - t[0].dotProduct(n);
		Vector3D projection = p.add(n.scalarMultiply(-dist));
		//endRegion

		Vector3D ap = projection.subtract(t[0]);

		//region nearest point is corners
		final double abDOTap = ab.dotProduct(ap);
		final double acDOTap = ac.dotProduct(ap);

		if (abDOTap <= 0d && acDOTap <= 0d) return t[0];

		final Vector3D bc = t[2].subtract(t[1]);
		final Vector3D bp = projection.subtract(t[1]);

		final double baDOTbp = ab.negate().dotProduct(bp);
		final double bcDOTbp = bc.dotProduct(bp);
		if (baDOTbp <= 0d && bcDOTbp <= 0d) return t[1];


		final Vector3D cp = projection.subtract(t[2]);
		final double cbDOTcp = bc.negate().dotProduct(cp);
		final double caDOTcp = ac.negate().dotProduct(cp);
		if (cbDOTcp <= 0d && caDOTcp <= 0d) return t[2];
		//endregion

		double acDOTac = ac.dotProduct(ac);
		double abDOTac =  ab.dotProduct(ac);
		double abDOTab = ab.dotProduct(ab);

		// Compute barycentric coordinates (v, w) of projection point
		double denom = (acDOTac * abDOTab - abDOTac *abDOTac);
		if (Math.abs(denom) < 1.0e-30) {
			return new Vector3D(-100,-100,-100); // Triangle is degenerate
		}

		double w = (acDOTac * abDOTap - abDOTac * acDOTap)/denom; //coordinate towards b from a
		double v = (abDOTab * acDOTap - abDOTac * abDOTap)/denom; //coordinate towards c from a


		// Check barycentric coordinates
		if ((v >= 0) && (w >= 0) && (v + w <= 1)) {
			// Nearest orthogonal projection point is in triangle
			return projection;
		}

		//region nearest point is on edge
		if(w <= 0 && v > w){
			return t[0].add(ab.scalarMultiply(v));
		}

		if(v <= 0 && w > v){
			return t[0].add(ac.scalarMultiply(w));
		}

		if(v + w > 1){
			final double scalarValue = bcDOTbp/bc.getNormSq();
			return t[1].add(bc.scalarMultiply(scalarValue));
		}
		//endregion

		if (v <=0 && w <= 0){ //this should be redundant, but for some reason isn't
			return t[0];
		}

		return new Vector3D(-100,-100,-100);
	}
}
