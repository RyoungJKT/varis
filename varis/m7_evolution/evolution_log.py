"""Evolution Log — Public, auditable record of every change Varis makes to itself.

Every action Varis takes on itself is recorded: retraining events, tool proposals,
integration decisions, metric changes, model versions. Published on VarisDB.

Uses SQLAlchemy for persistent storage with a SQLite or PostgreSQL backend.

Priority: 1 (build first — this is the backbone of self-evolution transparency)
"""
import json
import logging
from datetime import datetime, timezone
from typing import Any

from sqlalchemy import Column, DateTime, Integer, String, Text, create_engine, desc
from sqlalchemy.orm import Session, declarative_base, sessionmaker

logger = logging.getLogger(__name__)

# Event types
EVENT_RETRAIN_START = "RETRAIN_START"
EVENT_RETRAIN_COMPLETE = "RETRAIN_COMPLETE"
EVENT_DEPLOY = "DEPLOY"
EVENT_REJECT = "REJECT"
EVENT_ROLLBACK = "ROLLBACK"
EVENT_AUTO_RETRAIN = "AUTO_RETRAIN"
EVENT_ERROR = "ERROR"
EVENT_TOOL_DISCOVERY = "TOOL_DISCOVERY"
EVENT_LLM_ASSESSMENT = "LLM_ASSESSMENT"
EVENT_TOOL_INTEGRATION = "TOOL_INTEGRATION"

_VALID_EVENT_TYPES = {
    EVENT_RETRAIN_START,
    EVENT_RETRAIN_COMPLETE,
    EVENT_DEPLOY,
    EVENT_REJECT,
    EVENT_ROLLBACK,
    EVENT_AUTO_RETRAIN,
    EVENT_ERROR,
    EVENT_TOOL_DISCOVERY,
    EVENT_LLM_ASSESSMENT,
    EVENT_TOOL_INTEGRATION,
}

Base = declarative_base()


class EvolutionLogEntry(Base):
    """SQLAlchemy ORM model for evolution_log table."""

    __tablename__ = "evolution_log"

    id = Column(Integer, primary_key=True, autoincrement=True)
    event_type = Column(String, nullable=False)
    model_version = Column(String, nullable=True)
    details = Column(Text, nullable=True)
    timestamp = Column(DateTime, default=lambda: datetime.now(timezone.utc), index=True)

    def to_dict(self) -> dict[str, Any]:
        """Convert the ORM entry to a plain dict.

        Returns:
            Dict with id, event_type, model_version, details (parsed JSON), timestamp (ISO).
        """
        parsed_details = None
        if self.details is not None:
            try:
                parsed_details = json.loads(self.details)
            except (json.JSONDecodeError, TypeError):
                parsed_details = self.details

        return {
            "id": self.id,
            "event_type": self.event_type,
            "model_version": self.model_version,
            "details": parsed_details,
            "timestamp": self.timestamp.isoformat() if self.timestamp else None,
        }


def init_evolution_log(database_url: str) -> sessionmaker:
    """Create engine, session factory, and ensure the evolution_log table exists.

    Args:
        database_url: SQLAlchemy database URL (e.g. 'sqlite:///evo.db').

    Returns:
        A sessionmaker factory — call it to get a new Session.
    """
    try:
        engine = create_engine(database_url)
        Base.metadata.create_all(engine)
        factory = sessionmaker(bind=engine)
        logger.info("Evolution log initialized with database: %s", database_url)
        return factory
    except Exception as e:
        logger.error("Failed to initialize evolution log database: %s", e)
        raise


def log_event(
    session_factory: sessionmaker,
    event_type: str,
    model_version: str | None = None,
    details: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Record an evolution event to the log.

    Args:
        session_factory: Callable that returns a SQLAlchemy Session.
        event_type: One of the EVENT_* constants.
        model_version: Model version string (e.g. 'v2026.03'), or None.
        details: Dict with event-specific fields; serialized as JSON.

    Returns:
        The complete log entry as a dict with id, event_type, model_version,
        details, and timestamp.
    """
    if event_type not in _VALID_EVENT_TYPES:
        logger.warning("Unknown event type '%s'; recording anyway.", event_type)

    details_json = None
    if details is not None:
        try:
            details_json = json.dumps(details)
        except (TypeError, ValueError) as e:
            logger.warning("Could not serialize details to JSON: %s", e)
            details_json = json.dumps({"_serialization_error": str(e)})

    entry = EvolutionLogEntry(
        event_type=event_type,
        model_version=model_version,
        details=details_json,
        timestamp=datetime.now(timezone.utc),
    )

    try:
        session: Session = session_factory()
        try:
            session.add(entry)
            session.commit()
            session.refresh(entry)
            result = entry.to_dict()
            return result
        finally:
            session.close()
    except Exception as e:
        logger.error("Failed to log evolution event '%s': %s", event_type, e)
        raise


def get_log(
    session_factory: sessionmaker,
    limit: int = 50,
    event_type: str | None = None,
) -> list[dict[str, Any]]:
    """Retrieve evolution log entries, optionally filtered by type.

    Args:
        session_factory: Callable that returns a SQLAlchemy Session.
        limit: Maximum number of entries to return. Defaults to 50.
        event_type: If provided, filter to only this event type.

    Returns:
        List of dicts ordered by timestamp descending. Each dict has:
        id, event_type, model_version, details (parsed from JSON), timestamp (ISO string).
    """
    try:
        session: Session = session_factory()
        try:
            query = session.query(EvolutionLogEntry)
            if event_type is not None:
                query = query.filter(EvolutionLogEntry.event_type == event_type)
            query = query.order_by(desc(EvolutionLogEntry.timestamp))
            query = query.limit(limit)
            entries = query.all()
            return [entry.to_dict() for entry in entries]
        finally:
            session.close()
    except Exception as e:
        logger.error("Failed to retrieve evolution log: %s", e)
        return []


def get_current_model_version(session_factory: sessionmaker) -> str | None:
    """Return the current deployed model version string.

    Finds the most recent DEPLOY event and returns its model_version.

    Args:
        session_factory: Callable that returns a SQLAlchemy Session.

    Returns:
        The model_version string from the most recent DEPLOY event,
        or None if no DEPLOY events exist.
    """
    try:
        session: Session = session_factory()
        try:
            entry = (
                session.query(EvolutionLogEntry)
                .filter(EvolutionLogEntry.event_type == EVENT_DEPLOY)
                .order_by(desc(EvolutionLogEntry.timestamp))
                .first()
            )
            if entry is not None:
                return entry.model_version
            return None
        finally:
            session.close()
    except Exception as e:
        logger.error("Failed to get current model version: %s", e)
        return None
