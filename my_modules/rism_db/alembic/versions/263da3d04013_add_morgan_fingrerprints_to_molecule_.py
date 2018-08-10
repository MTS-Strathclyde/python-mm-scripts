"""add morgan fingrerprints to molecule table

Revision ID: 263da3d04013
Revises: None
Create Date: 2014-02-24 12:15:39.408295

"""

# revision identifiers, used by Alembic.
revision = '263da3d04013'
down_revision = None

from alembic import op
import sqlalchemy as sa

from razi.orm import ChemColumn
from razi.chemtypes import BitFingerprint

def upgrade():
    op.add_column('Molecule', ChemColumn('morgan', BitFingerprint))


def downgrade():
    op.drop_column('Molecule', 'morgan')
