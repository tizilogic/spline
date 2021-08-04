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
}
