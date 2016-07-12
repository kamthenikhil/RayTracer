/**
 * ray_tracer.cpp
 * CS230
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include <math.h>

#include "../Headers/ray_tracer.h"

using namespace std;

const double Object::small_t = 1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x) {
	return x * x;
}

Pixel Pixel_Color(const Vector_3D<double>& color) {
	Pixel pixel = 0;
	SET_RED(pixel, (unsigned char )(min(color.x, 1.0) * 255));
	SET_GREEN(pixel, (unsigned char )(min(color.y, 1.0) * 255));
	SET_BLUE(pixel, (unsigned char )(min(color.z, 1.0) * 255));
	return pixel;
}

//--------------------------------------------------------------------------------
// KD-Tree
//--------------------------------------------------------------------------------
KdTreeNode* KdTree::build(vector<Object*> objects) {
	/**
	 * Create a new node in the KD-Tree
	 */
	KdTreeNode* node = new KdTreeNode(objects);
	node->setMinimumBoundingBox();
	vector<Object*> objectsLeftNode;
	vector<Object*> objectsRightNode;
	if (objects.size() == 0) {
		return node;
	}
	if (objects.size() == 1) {
		node->leftChild = new KdTreeNode(objectsLeftNode);
		node->rightChild = new KdTreeNode(objectsRightNode);
		return node;
	}
	for (unsigned int i = 0; i < objects.size(); i++) {
		switch (node->minimumBoundingBox.longestAxis) {
		case 0:
			if (objects[i]->minimumBoundingBox.x_max
					> node->minimumBoundingBox.midpointValue) {
				objectsRightNode.push_back(objects[i]);
			} else {
				objectsLeftNode.push_back(objects[i]);
			}
			break;
		case 1:
			if (objects[i]->minimumBoundingBox.y_max
					> node->minimumBoundingBox.midpointValue) {
				objectsRightNode.push_back(objects[i]);
			} else {
				objectsLeftNode.push_back(objects[i]);
			}
			break;
		case 2:
			if (objects[i]->minimumBoundingBox.z_max
					> node->minimumBoundingBox.midpointValue) {
				objectsRightNode.push_back(objects[i]);
			} else {
				objectsLeftNode.push_back(objects[i]);
			}
			break;
		default:
			break;
		}
	}
	node->leftChild = build(objectsLeftNode);
	node->rightChild = build(objectsRightNode);
	return node;
}

void KdTreeNode::setMinimumBoundingBox() {
	minimumBoundingBox.x_max = -DOUBLE_MAX;
	minimumBoundingBox.y_max = -DOUBLE_MAX;
	minimumBoundingBox.z_max = -DOUBLE_MAX;
	minimumBoundingBox.x_min = DOUBLE_MAX;
	minimumBoundingBox.y_min = DOUBLE_MAX;
	minimumBoundingBox.z_min = DOUBLE_MAX;
	double cumulativeXMax = 0;
	double cumulativeYMax = 0;
	double cumulativeZMax = 0;
	long objectCount = objects.size();
	for(unsigned int i=0; i<objectCount; i++) {
		if(objects[i]->minimumBoundingBox.x_max > minimumBoundingBox.x_max) {
			minimumBoundingBox.x_max = objects[i]->minimumBoundingBox.x_max;
		}
		if(objects[i]->minimumBoundingBox.y_max > minimumBoundingBox.y_max) {
			minimumBoundingBox.y_max = objects[i]->minimumBoundingBox.y_max;
		}
		if(objects[i]->minimumBoundingBox.z_max > minimumBoundingBox.z_max) {
			minimumBoundingBox.z_max = objects[i]->minimumBoundingBox.z_max;
		}
		if(objects[i]->minimumBoundingBox.x_min < minimumBoundingBox.x_min) {
			minimumBoundingBox.x_min = objects[i]->minimumBoundingBox.x_min;
		}
		if(objects[i]->minimumBoundingBox.y_min < minimumBoundingBox.y_min) {
			minimumBoundingBox.y_min = objects[i]->minimumBoundingBox.y_min;
		}
		if(objects[i]->minimumBoundingBox.z_min < minimumBoundingBox.z_min) {
			minimumBoundingBox.z_min = objects[i]->minimumBoundingBox.z_min;
		}
		cumulativeXMax += objects[i]->minimumBoundingBox.x_max;
		cumulativeYMax += objects[i]->minimumBoundingBox.y_max;
		cumulativeZMax += objects[i]->minimumBoundingBox.z_max;
	}
	double lengthAlongX = minimumBoundingBox.x_max - minimumBoundingBox.x_min;
	double lengthAlongY = minimumBoundingBox.y_max - minimumBoundingBox.y_min;
	double lengthAlongZ = minimumBoundingBox.z_max - minimumBoundingBox.z_min;
	if(lengthAlongY > lengthAlongZ) {
		if(lengthAlongX > lengthAlongY) {
			minimumBoundingBox.longestAxis = 0;
			minimumBoundingBox.midpointValue = cumulativeXMax/objectCount;
		} else {
			minimumBoundingBox.longestAxis = 1;
			minimumBoundingBox.midpointValue = cumulativeYMax/objectCount;
		}
	} else {
		if(lengthAlongX > lengthAlongZ) {
			minimumBoundingBox.longestAxis = 0;
			minimumBoundingBox.midpointValue = cumulativeXMax/objectCount;
		} else {
			minimumBoundingBox.longestAxis = 2;
			minimumBoundingBox.midpointValue = cumulativeZMax/objectCount;
		}
	}
}

//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::Shade_Surface(const Ray& ray,
		const Object& intersection_object,
		const Vector_3D<double>& intersection_point,
		const Vector_3D<double>& same_side_normal) const {
	Vector_3D<double> color = Get_Ambient_Diffuse_Specular_Color(
			intersection_point, same_side_normal, ray);
	return color;
}

Vector_3D<double> Reflective_Shader::Shade_Surface(const Ray& ray,
		const Object& intersection_object,
		const Vector_3D<double>& intersection_point,
		const Vector_3D<double>& same_side_normal) const {
	// Fetch ambient, diffuse and specular components of color from Phong Shader.
	Vector_3D<double> color = Phong_Shader::Shade_Surface(ray,
			intersection_object, intersection_point, same_side_normal);
	Vector_3D<double> reflectiveColor(0, 0, 0);
	if (reflectivity > 0 && ray.recursion_depth > 0) {
		Vector_3D<double> reflectedRayDirection =
				Shader::Get_Reflected_Vector(
						Vector_3D<double>(0, 0, 0) - ray.direction,
						same_side_normal).Normalized();
		Ray reflectedRay(intersection_point, reflectedRayDirection);
		reflectedRay.recursion_depth = ray.recursion_depth - 1;
		reflectiveColor = world.Cast_Ray(reflectedRay, ray) * reflectivity;
	}
	color = color + reflectiveColor;
	return color;
}

Vector_3D<double> Refractive_Shader::Shade_Surface(const Ray& ray,
		const Object& intersection_object,
		const Vector_3D<double>& intersection_point,
		const Vector_3D<double>& same_side_normal) const {
	// Fetch ambient, diffuse, specular, reflective components of color from Reflectiive Shader.
	Vector_3D<double> color = Reflective_Shader::Shade_Surface(ray,
			intersection_object, intersection_point, same_side_normal);

	double n1;
	double n2;
	Vector_3D<double> normal;
	// Outside medium is considered to be air (with refractive Index: 1).
	// Depending on the angle between view ray and the normal we can decide
	// if the ray is falling on the object or leaving the object.
	// n1 and n2 are set accordingly.
	if (Vector_3D<double>::Dot_Product(ray.direction, same_side_normal) < 0) {
		n1 = 1;
		n2 = refractiveIndex;
		normal = same_side_normal;
	} else {
		n1 = refractiveIndex;
		n2 = 1;
		normal = Vector_3D<double>(0, 0, 0) - same_side_normal;
	}
	double refractionFactor = 1 - reflectivity;
	Vector_3D<double> refractiveColor(0, 0, 0);
	if (refractionFactor > 0) {
		Vector_3D<double> refractedRayDirection = Shader::Get_Refracted_Vector(
				ray.direction, normal, n1, n2).Normalized();
		Ray refractedRay(intersection_point, refractedRayDirection);
		refractiveColor = world.Cast_Ray(refractedRay, ray) * refractionFactor;
	}
	color = color + refractiveColor;
	return color;
}

// The following method is used to compute the Intensity of light transmitted, using Fresnel's Law.
double Refractive_Shader::Get_Reflectance(const Vector_3D<double>& direction,
		const Vector_3D<double>& normal, const double n1, const double n2) {
	double ratio = n1 / n2;
	double cosAngleOfIncidence = -Vector_3D<double>::Dot_Product(normal,
			direction);
	double cosAngleOfTransmittance = sqrt(
			1
					- (ratio * ratio
							* (1 - cosAngleOfIncidence * cosAngleOfIncidence)));
	double energyReflectedPerpendicular = (n1 * cosAngleOfIncidence
			- n2 * cosAngleOfTransmittance)
			/ (n1 * cosAngleOfIncidence + n2 * cosAngleOfTransmittance);
	double energyReflectedParallel = (n2 * cosAngleOfIncidence
			- n1 * cosAngleOfTransmittance)
			/ (n2 * cosAngleOfIncidence + n1 * cosAngleOfTransmittance);
	return pow(energyReflectedParallel, 2)
			+ pow(energyReflectedPerpendicular, 2);
}

Vector_3D<double> Flat_Shader::Shade_Surface(const Ray& ray,
		const Object& intersection_object,
		const Vector_3D<double>& intersection_point,
		const Vector_3D<double>& same_side_normal) const {
	return color;
}

// The following method is used to determine the ambient, diffuse and specular components of the color.
Vector_3D<double> Phong_Shader::Get_Ambient_Diffuse_Specular_Color(
		const Vector_3D<double>& intersection_point,
		const Vector_3D<double>& same_side_normal, const Ray& ray) const {
	Vector_3D<double> ambientColor(0.0, 0.0, 0.0);
	Vector_3D<double> diffuseColor(0.0, 0.0, 0.0);
	Vector_3D<double> specularColor(0.0, 0.0, 0.0);

	for (unsigned int i = 0; i < world.lights.size(); i++) {
		ambientColor = ambientColor
				+ color_ambient * world.lights[i]->color
						* world.lights[i]->brightness;
		Vector_3D<double> lightRay = intersection_point
				- world.lights[i]->position;
		Vector_3D<double> lightDirection = lightRay.Normalized();
		if (world.enable_shadows) {
			Ray shadowRay(intersection_point,
					Vector_3D<double>(0, 0, 0) - lightDirection);
			world.Closest_Intersection(shadowRay);
			if (!shadowRay.semi_infinite
					&& (lightRay.Length() > shadowRay.t_max)) {
				continue;
			}
		}
		double dotProductDiffuse = Vector_3D<double>::Dot_Product(
				Vector_3D<double>(0, 0, 0) - lightDirection, same_side_normal);
		Vector_3D<double> diffuseFactor(0, 0, 0);
		if (dotProductDiffuse > 0) {
			diffuseFactor = (color_diffuse * dotProductDiffuse)
					* (world.lights[i]->color) * world.lights[i]->brightness;
		}
		diffuseColor = diffuseColor + diffuseFactor;
		Vector_3D<double> view =
				(ray.endpoint - intersection_point).Normalized();

		Vector_3D<double> reflectedRay = Shader::Get_Reflected_Vector(
				Vector_3D<double>(0, 0, 0) - lightDirection, same_side_normal);

		double dotProductSpecular = Vector_3D<double>::Dot_Product(view,
				reflectedRay);
		Vector_3D<double> specularFactor(0, 0, 0);
		if (dotProductSpecular > 0) {
			specularFactor = color_specular
					* pow(dotProductSpecular, specular_power)
					* (world.lights[i]->color) * world.lights[i]->brightness;
		}
		specularColor = specularColor + specularFactor;
	}
	return ambientColor + diffuseColor + specularColor;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::Intersection(Ray& ray) const {
	double a = 1;
	double b = 2
			* Vector_3D<double>::Dot_Product(ray.direction,
					ray.endpoint - center);
	double c = Vector_3D<double>::Dot_Product(ray.endpoint - center,
			ray.endpoint - center) - pow(radius, 2);

	double determinant = pow(b, 2) - 4 * a * c;
	// Check if the the determinant is greater than 0, which implies
	if (determinant > small_t) {
		double t = 0;
		double t1 = (-b + sqrt(determinant)) / (2 * a);
		double t2 = (-b - sqrt(determinant)) / (2 * a);
		// As t1 > t2 the following checks will suffice
		if (t2 > small_t) {
			t = t2;
		} else {
			if (t1 > small_t) {
				t = t1;
			}
		}
		if (t > small_t) {
			ray.t_max = t;
			return true;
		}
	}
	return false;
}

Vector_3D<double> Sphere::Normal(const Vector_3D<double>& location) const {
	Vector_3D<double> normal((location - center).Normalized());
	return normal;
}

Vector_3D<double> Sphere::getMinimumBounds() const {
	Vector_3D<double> minimumBounds(center.x - radius, center.y - radius,
			center.z - radius);
	return minimumBounds;
}

Vector_3D<double> Sphere::getMaximumBounds() const {
	Vector_3D<double> maximumBounds(center.x + radius, center.y + radius,
			center.z + radius);
	return maximumBounds;
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::Intersection(Ray& ray) const {
	double denominator = Vector_3D<double>::Dot_Product(ray.direction, normal);
	if (abs(denominator) > small_t) {
		double t = Vector_3D<double>::Dot_Product(x1 - ray.endpoint, normal)
				/ denominator;
		if (t > small_t) {
			// Current_object and semi_infinite are set in Closest_Intersection().
			ray.t_max = t;
			return true;
		}
	}
	return false;
}

Vector_3D<double> Plane::Normal(const Vector_3D<double>& location) const {
	return normal;
}

Vector_3D<double> Plane::getMinimumBounds() const {
	Vector_3D<double> minimumBounds(0, 0, 0);
	return minimumBounds;
}

Vector_3D<double> Plane::getMaximumBounds() const {
	Vector_3D<double> maximumBounds(0, 0, 0);
	return maximumBounds;
}
//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel.
vector<Vector_3D<double> > Camera::World_Position(
		const Vector_2D<int>& pixel_index) {
	// Implemented Anti-Aliasing by considering 5 points per pixel. Each pixel is divided in 4 smaller squares.
	// Center of the pixel and center of the 4 smaller squares has been used to cast 5 view rays per pixel.
	vector<Vector_3D<double> > result;
	Vector_2D<double> center = film.pixel_grid.X(pixel_index);

	double dx = film.pixel_grid.dx;
	double dy = film.pixel_grid.dy;

	Vector_2D<double> leftTop(center.x - (dx / 4), center.y + (dy / 4));
	Vector_2D<double> leftBottom(center.x - (dx / 4), center.y - (dy / 4));
	Vector_2D<double> rightTop(center.x + (dx / 4), center.y + (dy / 4));
	Vector_2D<double> rightBottom(center.x + (dx / 4), center.y - (dy / 4));
//
	Vector_3D<double> centerPosition = focal_point
			+ horizontal_vector * center.x + vertical_vector * center.y;
	Vector_3D<double> leftTopPosition = focal_point
			+ horizontal_vector * leftTop.x + vertical_vector * leftTop.y;
	Vector_3D<double> leftBottomPosition = focal_point
			+ horizontal_vector * leftBottom.x + vertical_vector * leftBottom.y;
	Vector_3D<double> rightTopPosition = focal_point
			+ horizontal_vector * rightTop.x + vertical_vector * rightTop.y;
	Vector_3D<double> rightBottomPosition = focal_point
			+ horizontal_vector * rightBottom.x
			+ vertical_vector * rightBottom.y;

	result.push_back(centerPosition);
	result.push_back(leftTopPosition);
	result.push_back(leftBottomPosition);
	result.push_back(rightTopPosition);
	result.push_back(rightBottomPosition);

	return result;
}

bool Render_World::intersection(MinimumBoundingBox minimumBoundingBox,
		Ray& ray) {

	Vector_3D<double> directionFraction(1.0f / ray.direction.x,
			1.0f / ray.direction.y, 1.0f / ray.direction.z);

	double t1 = (minimumBoundingBox.x_min - ray.endpoint.x)
			* directionFraction.x;
	double t2 = (minimumBoundingBox.x_max - ray.endpoint.x)
			* directionFraction.x;
	double t3 = (minimumBoundingBox.y_min - ray.endpoint.y)
			* directionFraction.y;
	double t4 = (minimumBoundingBox.y_max - ray.endpoint.y)
			* directionFraction.y;
	double t5 = (minimumBoundingBox.z_min - ray.endpoint.z)
			* directionFraction.z;
	double t6 = (minimumBoundingBox.z_max - ray.endpoint.z)
			* directionFraction.z;

	double tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
	double tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

	/**
	 * This is the case when ray-MBR intersection happens in negative direction.
	 */
	if (tmax < 0) {
		return false;
	}
	/**
	 * This is the case when the ray doesn't intersect with the MBR.
	 */
	if (tmin > tmax) {
		return false;
	}
	return true;
}

bool Render_World::intersectionWithKdTreeNode(Ray& ray, KdTreeNode* node,
		double& t_max) {

	bool found = false;
	if (intersection(node->minimumBoundingBox, ray)) {
		if (node->leftChild->objects.size() > 0
				|| node->rightChild->objects.size() > 0) {
			bool leftIntersected = intersectionWithKdTreeNode(ray,
					node->leftChild, t_max);
			bool rightIntersected = intersectionWithKdTreeNode(ray,
					node->rightChild, t_max);
			return leftIntersected || rightIntersected;
		} else {
			for (unsigned int i = 0; i < node->objects.size(); i++) {
				found |= node->objects[i]->Intersection(ray);
				if (found && ray.t_max < t_max) {
					t_max = ray.t_max;
					ray.semi_infinite = false;
					ray.current_object = node->objects[i];
				}
			}
			return found;
		}
	}
	return found;
}

//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::Closest_Intersection(Ray& ray) {

	double t_max = FLT_MAX;
	if(useIndex) {
		if (intersectionWithKdTreeNode(ray, kdtree->root, t_max)) {
			ray.t_max = t_max;
		}
	} else {
		bool found = false;
		for (unsigned int i = 0; i < objects.size(); i++) {
			found |= objects[i]->Intersection(ray);
			if (found && ray.t_max < t_max) {
				t_max = ray.t_max;
				ray.current_object = objects[i];
			}
		}
		if (found) {
			ray.t_max = t_max;
			ray.semi_infinite = false;
			return ray.current_object;
		}
	}
	return 0;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const Vector_2D<int>& pixel_index) {
	Vector_3D<double> color(0, 0, 0);
	vector<Vector_3D<double> > results = camera.World_Position(pixel_index);
	double counter = 0;
	for (unsigned int i = 0; i < results.size(); i++) {
		Ray ray(camera.position, (results[i] - camera.position).Normalized());
		ray.recursion_depth = recursion_depth_limit;
		Ray dummy_root;
		color += Cast_Ray(ray, dummy_root);
		counter++;
	}
	Vector_3D<double> temp((color.x) / counter, (color.y) / counter,
			(color.z) / counter);
	camera.film.Set_Pixel(pixel_index, Pixel_Color(temp));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
Vector_3D<double> Render_World::Cast_Ray(Ray& ray, const Ray& parent_ray) {
	Flat_Shader* shader = dynamic_cast<Flat_Shader*>(background_shader);
	Vector_3D<double> color(shader->color);
	Closest_Intersection(ray);
	if (!ray.semi_infinite) {
		color = ray.current_object->material_shader->Shade_Surface(ray,
				*(ray.current_object), ray.Point(ray.t_max),
				ray.current_object->Normal(ray.Point(ray.t_max)));
	}
	return color;
}
