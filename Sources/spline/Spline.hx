package spline;

import kha.math.FastVector2;
import kha.FastFloat;


typedef Point2 = {x:FastFloat, y:FastFloat};
typedef OrientedPoint = {pos:Point2, fw:FastVector2};


class Spline {
    /**
     * Return a point on a Catmull Rom spline using 4 control points.
     * @param t Value in the range 0..1
     * @param p0 Control point for the start of the spline
     * @param p1 Start point of the spline
     * @param p2 End point of the spline
     * @param p3 Control point for the end of the spline
     * @return Point2
     */
    public static inline function splinePoint(t:FastFloat, p0:Point2, p1:Point2, p2:Point2, p3:Point2):Point2 {
        var tt:FastFloat = t * t;
		var ttt:FastFloat = tt * t;

		var q1:FastFloat = -ttt + 2.0 * tt - t;
		var q2:FastFloat = 3.0 * ttt - 5.0 * tt + 2.0;
		var q3:FastFloat = -3.0 * ttt + 4.0 * tt + t;
		var q4:FastFloat = ttt - tt;

		var tx:FastFloat = 0.5 * (p0.x * q1 + p1.x * q2 + p2.x * q3 + p3.x * q4);
		var ty:FastFloat = 0.5 * (p0.y * q1 + p1.y * q2 + p2.y * q3 + p3.y * q4);

		return {x:tx, y:ty};
    }

    /**
     * Return a oriented point on a Bezier Curve
     * @param t Value in the range 0..1
     * @param p0 Start point of the curve
     * @param p1 Control point for the start of the curve
     * @param p2 Control point for the end of the curve
     * @param p3 End point of the curve
     * @return OrientedPoint
     */
    public static inline function bezierOP(t:FastFloat, p0:Point2, p1:Point2, p2:Point2, p3:Point2):OrientedPoint {
        var a:Point2 = lerpPoint(t, p0, p1);
        var b:Point2 = lerpPoint(t, p1, p2);
        var c:Point2 = lerpPoint(t, p2, p3);

        var d:Point2 = lerpPoint(t, a, b);
        var e:Point2 = lerpPoint(t, b, c);
        return {pos:lerpPoint(t, d, e), fw:(new FastVector2(e.x - d.x, e.y - d.y)).normalized()};
    }

    static inline function lerpPoint(t:FastFloat, a:Point2, b:Point2):Point2 {
        return {
            x:a.x + (b.x - a.x) * t,
            y:a.y + (b.y - a.y) * t
        };
    }

    public static inline function cubicHermitePoint(t:FastFloat, points:Array<Point2>, tangents:Array<Point2>, ?knots:Array<FastFloat>):Point2 {
        if (t > 0.999999) return points[points.length - 1];
        var n = points.length;    // number or points / tangents / knots
        var i = 0;
        var i0 = -1;
        var i1 = -1;
        var scale:FastFloat = -1.0;

        if(knots != null) {
            // find knot interval for t
            while (i < n - 1) {
                if(t >= knots[i] && t <= knots[i+1]) {
                    break;
                }
            }

            if (i == n - 1) throw "Out of bounds";

            i0 = i;
            i1 = i + 1;
            var k0 = knots[i0];
            var k1 = knots[i1];
            scale = k1 - k0;

            t = (t - k0) / scale;
        }
        else {
            t = t * (n - 1); // rescale t to [0, n-1]
            i0 = Math.floor(t);
            i1 = i0 + 1;

            if (i0 > n - 1) throw "Out of bounds";
            if (i0 == n - 1) i1 = i0;

            scale = i1 - i0;

            t = (t - i0) / scale;
        }

        /*derivative:
          var t2 = t * t;
          var h00 = 6 * t2 - 6 * t;
          var h10 = 3 * t2 - 4 * t + 1;
          var h01 = - 6 * t2 + 6 * t;
          var h11 = 3 * t2 - 2 * t;
        */
        var t2 = t * t;
        var it = 1.0 - t;
        var it2 = it * it;
        var tt = 2.0 * t;
        var h00 = (1.0 + tt) * it2;
        var h10 = t * it2;
        var h01 = t2 * (3.0 - tt);
        var h11 = t2 * (t - 1.0);

        var tn0 = new FastVector2(tangents[i0].x - points[i0].x, tangents[i0].y - points[i0].y);
        var tn1 = new FastVector2(tangents[i1].x - points[i1].x, tangents[i1].y - points[i1].y);
        return {
            x: h00 * points[i0].x + h10 * tn0.x * scale + h01 * points[i1].x + h11 * tn1.x * scale,
            y: h00 * points[i0].y + h10 * tn0.y * scale + h01 * points[i1].y + h11 * tn1.y * scale
        };
    }
}
