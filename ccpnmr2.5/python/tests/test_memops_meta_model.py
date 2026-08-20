"""Smoke tests for memops.metamodel — the data-model backbone."""
import pytest

from memops.metamodel.ImpConstants import instance_level
from memops.metamodel.MetaModel import (
    MetaAttribute,
    MetaClass,
    MetaDataType,
    MetaPackage,
)


class TestImports:
    """Verify all major MetaModel classes import without error."""

    def test_import_meta_model_classes(self):
        assert MetaAttribute is not None
        assert MetaClass is not None
        assert MetaDataType is not None
        assert MetaPackage is not None

    def test_imp_constants(self):
        from memops.general.Constants import infinity

        assert isinstance(instance_level, str)
        assert len(instance_level) > 0
        assert infinity == -1


class TestAllowedTags:
    """Verify the allowedTags mapping from TaggedValues module."""

    def test_allowed_tags_structure(self):
        from memops.metamodel.TaggedValues import allowedTags

        assert isinstance(allowedTags, dict)
        assert "MetaModelElement" in allowedTags
