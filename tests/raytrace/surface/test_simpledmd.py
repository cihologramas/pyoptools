"""Minimal test suite for SimpleDMD surface."""

from math import isclose, radians, sqrt

import pytest
from pyoptools.raytrace.ray.ray import Ray
from pyoptools.raytrace.shape.circular import Circular
from pyoptools.raytrace.shape.rectangular import Rectangular
from pyoptools.raytrace.surface.simpledmd import SimpleDMD

# ============================================================================
# Section 1: Instantiation and Basic Validation
# ============================================================================


def test_simpledmd_instantiation():
    """Test SimpleDMD can be instantiated with all parameters."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(45),
        off_direction_angle=radians(225),
        state="on",
        shape=Rectangular(size=(10, 10)),
    )
    assert dmd.state == "on"
    assert dmd.shape is not None


def test_simpledmd_default_state():
    """Test SimpleDMD defaults to 'flat' state."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(0),
        off_direction_angle=radians(180),
    )
    assert dmd.state == "flat"


def test_simpledmd_invalid_state():
    """Test SimpleDMD rejects invalid state values."""
    with pytest.raises(ValueError):
        SimpleDMD(
            tilt_angle=radians(12),
            on_direction_angle=radians(0),
            off_direction_angle=radians(180),
            state="invalid",
        )


# ============================================================================
# Section 2: State Transitions
# ============================================================================


def test_simpledmd_state_transitions():
    """Test state can be changed between valid values."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(0),
        off_direction_angle=radians(180),
        state="flat",
    )

    dmd.state = "on"
    assert dmd.state == "on"

    dmd.state = "off"
    assert dmd.state == "off"

    dmd.state = "flat"
    assert dmd.state == "flat"


def test_simpledmd_invalid_state_transition():
    """Test setting invalid state raises ValueError."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(0),
        off_direction_angle=radians(180),
    )

    with pytest.raises(ValueError):
        dmd.state = "broken"


# ============================================================================
# Section 3: Normal Vector Calculations
# ============================================================================


def test_simpledmd_cardinal_directions_normals():
    """
    For tilt=30°, test (on/off) direction angles of 0, 90, 180, 270 deg:
    - At each position, normal should match precomputed, fixed expected values (regression, not math test)
    """
    from math import isclose, radians

    tilt_deg = 30
    tilt = radians(tilt_deg)
    angles_deg = [0, 90, 180, 270]
    angles_rad = [radians(deg) for deg in angles_deg]

    # Precomputed normals for 30 degree tilt at four cardinal directions
    expected_normals = [
        (0.5, 0.0, 0.86602540378),  # 0 degrees
        (0.0, 0.5, 0.86602540378),  # 90 degrees
        (-0.5, 0.0, 0.86602540378),  # 180 degrees
        (0.0, -0.5, 0.86602540378),  # 270 degrees
    ]

    # ON state: check returned normal matches fixed values
    for i, on_dir in enumerate(angles_rad):
        dmd = SimpleDMD(
            tilt_angle=tilt,
            on_direction_angle=on_dir,
            off_direction_angle=angles_rad[(i + 2) % 4],  # Opposite dir
            state="on",
        )
        normal = dmd.normal((0, 0, 0))
        expected = expected_normals[i]
        for j, comp in enumerate(["x", "y", "z"]):
            assert isclose(normal[j], expected[j], abs_tol=1e-8), (
                f"ON dir angle {angles_deg[i]}°: {comp} = {normal[j]} (expected {expected[j]})"
            )
        # Only one component (x or y) should be nonzero (or both approximately zero)
        assert isclose(normal[0] * normal[1], 0.0, abs_tol=1e-9)

    # OFF state: same normals but for off_direction_angle
    for i, off_dir in enumerate(angles_rad):
        dmd = SimpleDMD(
            tilt_angle=tilt,
            on_direction_angle=angles_rad[(i + 2) % 4],
            off_direction_angle=off_dir,
            state="off",
        )
        normal = dmd.normal((0, 0, 0))
        expected = expected_normals[i]
        for j, comp in enumerate(["x", "y", "z"]):
            assert isclose(normal[j], expected[j], abs_tol=1e-8), (
                f"OFF dir angle {angles_deg[i]}°: {comp} = {normal[j]} (expected {expected[j]})"
            )
        assert isclose(normal[0] * normal[1], 0.0, abs_tol=1e-9)


def test_simpledmd_flat_state_normal():
    """Test flat state produces perpendicular normal (0, 0, 1)."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(45),
        off_direction_angle=radians(225),
        state="flat",
    )

    # Get normal at origin
    normal = dmd.normal((0, 0, 0))

    # Should be perpendicular regardless of tilt_angle
    assert isclose(normal[0], 0.0, abs_tol=1e-9)
    assert isclose(normal[1], 0.0, abs_tol=1e-9)
    assert isclose(normal[2], 1.0, abs_tol=1e-9)


def test_simpledmd_on_state_normal():
    """Test on state produces correctly tilted normal."""
    tilt = radians(12)
    direction = radians(0)  # Along +X axis

    dmd = SimpleDMD(
        tilt_angle=tilt,
        on_direction_angle=direction,
        off_direction_angle=radians(180),
        state="on",
    )

    normal = dmd.normal((0, 0, 0))

    # Expected normal for tilt along +X axis
    from math import cos, sin

    expected_x = sin(tilt) * cos(direction)
    expected_y = sin(tilt) * sin(direction)
    expected_z = cos(tilt)

    assert isclose(normal[0], expected_x, abs_tol=1e-9)
    assert isclose(normal[1], expected_y, abs_tol=1e-9)
    assert isclose(normal[2], expected_z, abs_tol=1e-9)

    # Verify unit length
    length = sqrt(normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2)
    assert isclose(length, 1.0, abs_tol=1e-9)


def test_simpledmd_off_state_normal():
    """Test off state produces correctly tilted normal in opposite direction."""
    tilt = radians(12)
    on_dir = radians(45)
    off_dir = radians(225)  # Opposite of 45°

    dmd = SimpleDMD(
        tilt_angle=tilt,
        on_direction_angle=on_dir,
        off_direction_angle=off_dir,
        state="off",
    )

    normal = dmd.normal((0, 0, 0))

    # Expected normal
    from math import cos, sin

    expected_x = sin(tilt) * cos(off_dir)
    expected_y = sin(tilt) * sin(off_dir)
    expected_z = cos(tilt)

    assert isclose(normal[0], expected_x, abs_tol=1e-9)
    assert isclose(normal[1], expected_y, abs_tol=1e-9)
    assert isclose(normal[2], expected_z, abs_tol=1e-9)

    # Verify unit length
    length = sqrt(normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2)
    assert isclose(length, 1.0, abs_tol=1e-9)


# ============================================================================
# Section 4: Ray Reflection (Basic Verification)
# ============================================================================


def test_simpledmd_reflects_ray_in_flat_state():
    """Test that a ray reflects off SimpleDMD in flat state (normal incidence)."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(0),
        off_direction_angle=radians(180),
        state="flat",
        shape=Circular(radius=10),
    )

    # Ray traveling in -Z direction (toward surface at Z=0)
    ray = Ray(origin=(0, 0, 5), direction=(0, 0, -1))

    # Find intersection
    intersection = dmd.intersection(ray)

    # Should intersect at Z=0
    assert isclose(intersection[0], 0.0, abs_tol=1e-6)
    assert isclose(intersection[1], 0.0, abs_tol=1e-6)
    assert isclose(intersection[2], 0.0, abs_tol=1e-6)

    # Propagate ray to surface
    # Get reflected ray (passing a plausible refractive index - only reflection should occur)
    reflected_rays = dmd.propagate(ray, 1.0, 1.0)  # Reflection only

    # Should have exactly one reflected ray
    assert len(reflected_rays) == 1


def test_simpledmd_normal_changes_with_state():
    """Test that normal vector changes when state changes."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(0),
        off_direction_angle=radians(180),
        state="flat",
    )

    normal_flat = dmd.normal((0, 0, 0))

    dmd.state = "on"
    normal_on = dmd.normal((0, 0, 0))

    dmd.state = "off"
    normal_off = dmd.normal((0, 0, 0))

    # All three should be different
    assert not (
        isclose(normal_flat[0], normal_on[0], abs_tol=1e-9)
        and isclose(normal_flat[1], normal_on[1], abs_tol=1e-9)
        and isclose(normal_flat[2], normal_on[2], abs_tol=1e-9)
    )

    assert not (
        isclose(normal_on[0], normal_off[0], abs_tol=1e-9)
        and isclose(normal_on[1], normal_off[1], abs_tol=1e-9)
        and isclose(normal_on[2], normal_off[2], abs_tol=1e-9)
    )


# ============================================================================
# Section 5: Representation
# ============================================================================


def test_simpledmd_repr():
    """Test SimpleDMD __repr__ returns expected format."""
    dmd = SimpleDMD(
        tilt_angle=radians(12),
        on_direction_angle=radians(45),
        off_direction_angle=radians(225),
        state="on",
        shape=Rectangular(size=(10, 10)),
    )

    repr_str = repr(dmd)

    # Should contain class name and key parameters
    assert "SimpleDMD" in repr_str
    assert "state='on'" in repr_str
    assert "tilt_angle" in repr_str
    assert "on_direction_angle" in repr_str
    assert "off_direction_angle" in repr_str
